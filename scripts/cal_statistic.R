library(dplyr)
source('rtest1grp.R')

sc_overlapped_path <- '../data/sc_adata_high_variable_genes.csv'
sc_labels_path <- '../data/sc_labels.csv'
output_path <- '../data/marker_genes.csv'

sc_overlapped <- read.csv(sc_overlapped_path, header = TRUE, row.names = 1, check.names = FALSE)
sc_anno <- read.csv(sc_labels_path, header = TRUE, sep = ',', row.names = 1, check.names = FALSE)

sc_anno <- sc_anno %>% arrange(cell_type)
sc_overlapped <- sc_overlapped[rownames(sc_anno), ]

df_t_test <- nrow(sc_overlapped) - 2 
df_chisqstat_test <- length(unique(sc_anno$cell_type)) - 2 


robust_quasi <- function(y, group, max2 = FALSE, ...) {
  y <- unlist(y)
  if (max2) {
    run <- rtest2grp(y, group, ...)
  } else {
    run <- rtest(y, group, ...)
  }
  
  if (run$rcode == 1) {
    return(c(run$out, run$b))
  } else {
    return(c(rep(NA, ifelse(max2, 5, 3)), run$b))
  }
}

# main
all_results <- list()

for (col_name in names(sc_overlapped)) {
  print(paste("Processing column:", col_name))
  result <- tryCatch({
    x <- split(sc_overlapped[[col_name]], as.factor(sc_anno$cell_type))
    m1 <- sapply(x, median)
    pick <- which(m1 > 0)
    p1 <- sc_anno$cell_type %in% names(m1[pick])
    gr <- sc_anno$cell_type[p1]
    x1 <- sc_overlapped[[col_name]][p1]
    if(length(unique(gr)) > 1 && any(gr > 2)){
      robust_quasi(x1, as.factor(gr), max2 = FALSE)[1:3]
    } else {
      robust_quasi(sc_overlapped[[col_name]], as.factor(sc_anno$cell_type), max2 = FALSE, espXWX=1e-8)[1:3]
    }
  }, error = function(e) {
    print(paste("Error in column:", col_name, " - ", e$message))
    return(NULL)
  })
  all_results[[col_name]] <- result
}

result_df <- do.call(rbind, all_results)
colnames(result_df) <- c("maxgroup", "tstat", "chisqstat")
result_df <- as.data.frame(result_df)

result_df_without_NaN <- result_df[result_df$maxgroup %in% sc_anno$cell_type, ]
result_df_without_NaN$maxgroup <- as.factor(as.character(result_df_without_NaN$maxgroup))

result_df_without_NaN$tstat <- as.numeric(result_df_without_NaN$tstat)
result_df_without_NaN$chisqstat <- as.numeric(result_df_without_NaN$chisqstat)

# t-test
t_statistic <- result_df_without_NaN$tstat
p_value_t_test <- 2 * pt(-abs(t_statistic), df_t_test)
result_df_without_NaN$p_value_t_test <- p_value_t_test

# Chisquare test
chisqstat_statistic <- result_df_without_NaN$chisqstat
p_value_chisqstat_test <- pchisq(chisqstat_statistic, df_chisqstat_test, lower.tail = FALSE)
result_df_without_NaN$p_value_chisqstat_test <- p_value_chisqstat_test

result_df_without_NaN$Gene <- row.names(result_df_without_NaN)

# Selecting for significant marker genes
result_df_top10 <- result_df_without_NaN %>%
  group_by(maxgroup) %>%
  arrange(desc(tstat)) %>%
  slice_head(n = 10) %>%
  ungroup()

result_df_top10 <- as.data.frame(result_df_top10)
rownames(result_df_top10) <- result_df_top10$Gene

write.csv(result_df_top10, file = output_path, row.names = FALSE)
