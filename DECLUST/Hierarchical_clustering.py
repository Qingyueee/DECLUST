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

def selected_number_of_clusters(st_highly_variable_genes_df):
    st_np = st_highly_variable_genes_df.values
    inertia_values = []
    silhouette_scores = []
    cluster_range = range(2, 11)

    for n_clusters in cluster_range:
        labels, inertia = fit_and_predict(st_np, n_clusters)
        inertia_values.append(inertia)
        silhouette_scores.append(silhouette_score(st_np, labels))

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

def hierarchical_clustering(st_highly_variable_genes_df, x, coords):
    st_np = st_highly_variable_genes_df.values
    labels, inertia = fit_and_predict(st_np, x)

    hierarchical_label_df = pd.DataFrame(labels, columns=['label']).set_index(st_highly_variable_genes_df.index)
    hierarchical_results = pd.merge(hierarchical_label_df, coords, left_index=True, right_index=True)

    unique_labels = sorted(hierarchical_results['label'].unique())
    num_clusters = len(unique_labels)

    cmap = plt.cm.get_cmap('tab20', num_clusters)

    plt.figure(figsize=(8, 6))
    for i, label in enumerate(unique_labels):
        subset = hierarchical_results[hierarchical_results['label'] == label]
        color = cmap(i)
        plt.scatter(subset['x'], subset['y'], label=f"Cluster {label}", color=color, s=8)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='x-small')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Hierarchical Results')
    plt.tight_layout()
    plt.show()

    return hierarchical_results