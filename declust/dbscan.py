import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import DBSCAN

def plot_dbscan_clusters(cluster_arr, dbscan_labels, hierarchical_label):
    outlier_indices = np.where(dbscan_labels == -1)[0]
    non_outlier_indices = np.where(dbscan_labels != -1)[0]
    x_outliers, y_outliers = cluster_arr[outlier_indices, 0], cluster_arr[outlier_indices, 1]
    x_non_outliers, y_non_outliers = cluster_arr[non_outlier_indices, 0], cluster_arr[non_outlier_indices, 1]
    
    plt.figure(figsize=(8, 6), dpi=80)
    plt.scatter(x_outliers, y_outliers, color='red', marker='o', label='Outliers', s=10)

    unique_labels = np.unique(dbscan_labels)
    unique_labels = unique_labels[unique_labels != -1]
    colors = plt.cm.jet(np.linspace(0, 1, len(unique_labels)))

    for i, label in enumerate(unique_labels):
        cluster_indices = np.where(dbscan_labels == label)[0]
        plt.scatter(cluster_arr[cluster_indices, 0], cluster_arr[cluster_indices, 1], 
                    color=colors[i], marker='o', label=f'Cluster {label}', s=8)

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
    plt.title(f'DBSCAN Subclusters in Hierarchical Cluster {hierarchical_label}')
    plt.tight_layout()
    plt.show()

def find_cluster_centers(cluster_arr, dbscan_labels, hierarchical_label, min_cluster_size_ratio=0.05, sample_ratio=0.5):
    results = []
    unique_dbscan_labels = np.unique(dbscan_labels)
    
    for dbscan_label in unique_dbscan_labels:
        if dbscan_label != -1:
            cluster_indices = np.where(dbscan_labels == dbscan_label)[0]
            cluster_size = len(cluster_indices)

            if cluster_size < len(cluster_arr) * min_cluster_size_ratio:
                continue
            
            # Calculate initial seeds
            cluster_center = np.mean(cluster_arr[cluster_indices], axis=0)
            nearest_point_index = np.argmin(np.linalg.norm(cluster_arr[cluster_indices] - cluster_center, axis=1))
            nearest_point = cluster_arr[cluster_indices][nearest_point_index]
            
            points_to_add = [{'center_x': nearest_point[0], 'center_y': nearest_point[1]}]

            sample_size = max(1, int(cluster_size * sample_ratio))
            sample_indices = np.random.choice(cluster_indices, size=sample_size, replace=False)
            sampled_points = cluster_arr[sample_indices]
            
            for point in sampled_points:
                point_dict = {'center_x': point[0], 'center_y': point[1]}
                if point_dict not in points_to_add:
                    points_to_add.append(point_dict)
            
            for point_dict in points_to_add:
                results.append({
                    'hierarchical_label': hierarchical_label,
                    'dbscan_label': dbscan_label,
                    'center_x': point_dict['center_x'],
                    'center_y': point_dict['center_y']
                })
    
    return results

def clustering(hierarchical_results, coords, visualize=False):

    """
    Performs DBSCAN clustering on the results from hierarchical clustering and returns 
    the cluster centers mapped to the original coordinates indices.

    Parameters:
        hierarchical_results (pandas.DataFrame):
            DataFrame containing the hierarchical clustering results. It must include at least:
                - 'label': The label for each data point from the hierarchical clustering.
                - 'x', 'y': The coordinates of the data points.
                
        coords (pandas.DataFrame):
            DataFrame of the original coordinates. It should have 'x' and 'y' columns that 
            will be used to map the DBSCAN cluster centers back to the original indices.
            
        visualize (bool, optional):
            Whether to visualize the DBSCAN clustering for each cluster using the 
            `plot_dbscan_clusters()` function. Defaults to False.

    Returns:
        pandas.DataFrame:
            A DataFrame containing the DBSCAN cluster centers. The index of the DataFrame
            corresponds to the indices in the original `coords` DataFrame, and it includes
            center coordinates (e.g., 'center_x', 'center_y') for each cluster.
    """

    print("Running DBSCAN clustering...")
    eps = 4
    min_samples = 8
    results = []
    unique_cluster_labels = hierarchical_results['label'].unique()

    for cluster_label in unique_cluster_labels:
        cluster_data = hierarchical_results[hierarchical_results['label'] == cluster_label].copy()
        cluster_arr = cluster_data[['x', 'y']].values

        dbscan_labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(cluster_arr)

        if visualize:
            plot_dbscan_clusters(cluster_arr, dbscan_labels, cluster_label)
        
        centers = find_cluster_centers(cluster_arr, dbscan_labels, cluster_label)
        results.extend(centers)

    dbscan_centers_df = pd.DataFrame(results)
    coord_to_index = {(row.x, row.y): idx for idx, row in coords.iterrows()}
    matched_indexes = dbscan_centers_df.apply(
        lambda row: coord_to_index.get((row['center_x'], row['center_y']), None), axis=1
    )

    if matched_indexes.isnull().any():
        print("Warning: Some DBSCAN centers could not be matched to coords!")

    dbscan_centers_df.index = matched_indexes
    return dbscan_centers_df



