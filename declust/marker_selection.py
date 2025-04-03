from rpy2.robjects import r, pandas2ri
import importlib.resources
import pandas as pd

pandas2ri.activate()

def generate_marker_genes(sc_overlapped_path, sc_labels_path, output_path):

    """
    Generates marker genes by processing single-cell data through an R script using rpy2.
    
    This function reads two CSV files: one containing single-cell overlapped data (sc_overlapped)
    and one containing cell type annotations (sc_labels). It then calls an embedded R script to perform 
    differential analysis between cell types, computes t-test and chi-square statistics, and selects 
    the top 10 marker genes for each cell type. The results are saved to an output CSV file, which is 
    then read into a pandas DataFrame and returned.
    
    Parameters:
        sc_overlapped_path (str):
            Path to the CSV file containing the single-cell overlapped data. This file should have 
            row names corresponding to cell identifiers and include gene expression values.
            
        sc_labels_path (str):
            Path to the CSV file containing the cell type annotations. This file should have a column 
            named 'cell_type' and row names corresponding to cell identifiers.
            
        output_path (str):
            Path where the resulting CSV file with marker genes will be saved.
    
    Returns:
        pandas.DataFrame:
            A DataFrame containing the top marker genes per cell type, including columns for 
            maximum group, t-statistic, chi-square statistic, and associated p-values, with the 
            gene names as the index.
    """

    with importlib.resources.path("declust.scripts", "rtest1grp.R") as r_script_path:
        r_script_str = str(r_script_path)
        
    r(f'''
    library(dplyr)
    source('{r_script_str}')

    sc_overlapped <- read.csv('{sc_overlapped_path}', header = TRUE, row.names = 1, check.names = FALSE)
    sc_anno <- read.csv('{sc_labels_path}', header = TRUE, sep = ',', row.names = 1, check.names = FALSE)

    sc_anno <- sc_anno %>% arrange(cell_type)
    sc_overlapped <- sc_overlapped[rownames(sc_anno), ]

    df_t_test <- nrow(sc_overlapped) - 2 
    df_chisqstat_test <- length(unique(sc_anno$cell_type)) - 2 

    robust_quasi <- function(y, group, max2 = FALSE, ...) {{
      y <- unlist(y)
      if (max2) {{
        run <- rtest2grp(y, group, ...)
      }} else {{
        run <- rtest(y, group, ...)
      }}

      if (run$rcode == 1) {{
        return(c(run$out, run$b))
      }} else {{
        return(c(rep(NA, ifelse(max2, 5, 3)), run$b))
      }}
    }}

    all_results <- list()
    for (col_name in names(sc_overlapped)) {{
      result <- tryCatch({{
        x <- split(sc_overlapped[[col_name]], as.factor(sc_anno$cell_type))
        m1 <- sapply(x, median)
        pick <- which(m1 > 0)
        p1 <- sc_anno$cell_type %in% names(m1[pick])
        gr <- sc_anno$cell_type[p1]
        x1 <- sc_overlapped[[col_name]][p1]
        if(length(unique(gr)) > 1 && any(gr > 2)){{
          robust_quasi(x1, as.factor(gr), max2 = FALSE)[1:3]
        }} else {{
          robust_quasi(sc_overlapped[[col_name]], as.factor(sc_anno$cell_type), max2 = FALSE, espXWX=1e-8)[1:3]
        }}
      }}, error = function(e) {{
        return(NULL)
      }})
      all_results[[col_name]] <- result
    }}

    result_df <- do.call(rbind, all_results)
    colnames(result_df) <- c("maxgroup", "tstat", "chisqstat")
    result_df <- as.data.frame(result_df)

    result_df_without_NaN <- result_df[result_df$maxgroup %in% sc_anno$cell_type, ]
    result_df_without_NaN$maxgroup <- as.factor(as.character(result_df_without_NaN$maxgroup))

    result_df_without_NaN$tstat <- as.numeric(result_df_without_NaN$tstat)
    result_df_without_NaN$chisqstat <- as.numeric(result_df_without_NaN$chisqstat)

    result_df_without_NaN$chisqstat[result_df_without_NaN$chisqstat <= 0] <- 0.1
    result_df_without_NaN$log2chisqstat=log2(result_df_without_NaN$chisqstat+1)
    result_df_without_NaN$log2tstat=log2(result_df_without_NaN$tstat+1)

    t_statistic <- result_df_without_NaN$tstat
    p_value_t_test <- 2 * pt(-abs(t_statistic), df_t_test)
    result_df_without_NaN$p_value_t_test <- p_value_t_test
    result_df_without_NaN$FDR_t_test <- p.adjust(result_df_without_NaN$p_value_t_test, method = "BH")


    chisqstat_statistic <- result_df_without_NaN$chisqstat
    p_value_chisqstat_test <- pchisq(chisqstat_statistic, df_chisqstat_test, lower.tail = FALSE)
    result_df_without_NaN$p_value_chisqstat_test <- p_value_chisqstat_test
    result_df_without_NaN$FDR_chisqstat_test <- p.adjust(result_df_without_NaN$p_value_chisqstat_test, method = "BH")

    result_df_without_NaN$Gene <- row.names(result_df_without_NaN)

    pick=result_df_without_NaN$FDR_t_test < 0.05 & result_df_without_NaN$FDR_chisqstat_test > 0.05
    result_df_keep=result_df_without_NaN[pick,]
    # table(pick)
    # table(result_df_keep$maxgroup)
    # result_df_keep <- result_df_keep %>% arrange(desc(maxgroup))


    result_df_top10 <- result_df_without_NaN %>%
      group_by(maxgroup) %>%
      arrange(desc(tstat)) %>%
      slice_head(n = 10) %>%
      ungroup()

    result_df_top10 <- as.data.frame(result_df_top10)
    rownames(result_df_top10) <- result_df_top10$Gene
    combined_df <- rbind(result_df_top10, result_df_keep)
    marker_df <- unique(combined_df)
    marker_df <- marker_df %>% arrange(desc(maxgroup))

    write.csv(marker_df, file = '{output_path}', row.names = FALSE)
    ''')

    result_df = pd.read_csv(output_path, index_col='Gene')

    return result_df
