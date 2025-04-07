import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering


def fit_and_predict(X, n_clusters):
    model = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
    labels = model.fit_predict(X)
    inertia = calculate_inertia(X, labels)
    return labels, inertia


def calculate_inertia(X, labels):
    inertia = 0
    for label in np.unique(labels):
        cluster_mean = np.mean(X[labels == label], axis=0)
        inertia += np.sum((X[labels == label] - cluster_mean) ** 2)
    return inertia

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score

def selected_number_of_clusters(st_np, user_defined_clusters=None, show_selection_plot=False):

    inertia_values = []
    silhouette_scores = []
    cluster_range = range(2, 11)

    for n_clusters in cluster_range:
        labels, inertia = fit_and_predict(st_np, n_clusters)
        inertia_values.append(inertia)
        silhouette_scores.append(silhouette_score(st_np, labels))

    if show_selection_plot:
        fig, axs = plt.subplots(1, 2, figsize=(10, 4))
        axs[0].plot(cluster_range, inertia_values, marker='o')
        axs[0].set_title('Elbow Method for Optimal Number of Clusters')
        axs[0].set_xlabel('Number of clusters')
        axs[0].set_ylabel('Within-Cluster Sum of Squares (WCSS)')
        axs[0].set_xticks(cluster_range)

        axs[1].plot(cluster_range, silhouette_scores, marker='o', color='orange')
        axs[1].set_title('Silhouette Scores for Different Numbers of Clusters')
        axs[1].set_xlabel('Number of clusters')
        axs[1].set_ylabel('Silhouette Score')
        axs[1].set_xticks(cluster_range)

        plt.tight_layout()
        plt.show()

    inertia_diffs = np.diff(inertia_values)
    slope_diff_threshold = np.median(np.abs(inertia_diffs)) * 0.1 

    slowdown_index = None
    for i in range(1, len(inertia_diffs) - 1):
        slope_diff = abs(inertia_diffs[i] - inertia_diffs[i + 1])
        if slope_diff < slope_diff_threshold:
            slowdown_index = i + 1 
            break

    if slowdown_index is not None:
        silhouette_tail = silhouette_scores[slowdown_index:]
        if silhouette_tail:
            best_idx = np.argmax(silhouette_tail)
            auto_clusters = cluster_range[slowdown_index + best_idx]
        else:
            auto_clusters = cluster_range[np.argmax(silhouette_scores)]
    else:
        inertia_diffs2 = np.diff(inertia_values, 2)
        auto_clusters = cluster_range[np.argmin(inertia_diffs2) + 1]

    if user_defined_clusters is not None and user_defined_clusters > 1:
        return user_defined_clusters
    else:
        print(f"Using automatically selected number of clusters: {auto_clusters}")
        return auto_clusters


def clustering(st_highly_variable_genes_df, coords, num_clusters=None, show_plot=False, save_path=None, show_selection_plot=False):
    """
    Perform hierarchical clustering on spatial transcriptomics data based on highly variable genes.
    """
    import os

    print("Running Hierarchical Clustering...")
    st_np = st_highly_variable_genes_df.values

    n_clusters = selected_number_of_clusters(
        st_np,
        user_defined_clusters=num_clusters,
        show_selection_plot=show_selection_plot
    )

    labels, inertia = fit_and_predict(st_np, n_clusters)

    hierarchical_label_df = pd.DataFrame(labels, columns=['label']).set_index(st_highly_variable_genes_df.index)
    hierarchical_results = pd.merge(hierarchical_label_df, coords, left_index=True, right_index=True)

    unique_labels = sorted(hierarchical_results['label'].unique())
    cmap = plt.cm.get_cmap('tab20', len(unique_labels))

    plt.figure(figsize=(8, 6))
    for i, label in enumerate(unique_labels):
        subset = hierarchical_results[hierarchical_results['label'] == label]
        plt.scatter(subset['x'], subset['y'], label=f"Cluster {label}", color=cmap(i), s=8)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Hierarchical Clustering Result')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='x-small')
    plt.tight_layout()

    if show_plot:
        plt.show()
    elif save_path:
        os.makedirs(save_path, exist_ok=True)
        clustering_plot_file = os.path.join(save_path, "hierarchical_clustering_result.png")
        plt.savefig(clustering_plot_file, dpi=300)
        print(f"Cluster plot saved to: {clustering_plot_file}")
        plt.close()
    else:
        plt.close()

    return hierarchical_results
