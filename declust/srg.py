import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm
from sklearn.neighbors import KDTree
from sklearn.metrics.pairwise import euclidean_distances

def cluster_label(dbscan_centers_df):
    label_dict = {'index': [], 'coordinates': [], 'label': []}
    for idx, row in dbscan_centers_df.iterrows():
        label_dict['index'].append(idx)
        coordinates = (row['center_x'], row['center_y'])
        label_dict['coordinates'].append(coordinates)
        label_dict['label'].append(row['hierarchical_label'])
    return label_dict

def calculate_initial_distances(x, y, coords, gene_expression, directions, label_dict):

    result_dict = {'index': [], 'neighbour_coordinates': [], 'distances': [], 'initial_seed': []}

    initial_seed = (x, y)

    for direction in directions:
        distances = []

        dx, dy = direction
        neighbour_x, neighbour_y = x + dx, y + dy

        if (neighbour_x, neighbour_y) not in zip(coords['x'], coords['y']):
            continue

        if (neighbour_x, neighbour_y) in label_dict['coordinates']:
            continue

        neighbour_index = (np.abs(coords['x'] - neighbour_x) + np.abs(coords['y'] - neighbour_y)).idxmin()
        neighbour_gene_expression = gene_expression.loc[neighbour_index]

        distance = np.linalg.norm(gene_expression.loc[coords[(coords['x'] == x) & (coords['y'] == y)].index[0]] - neighbour_gene_expression)

        distances.append(distance)

        result_dict['index'].append(neighbour_index)
        result_dict['neighbour_coordinates'].append((x + dx, y + dy))
        result_dict['distances'].append(distances[0])  # distances only contains the result of one iteration
        result_dict['initial_seed'].append(initial_seed)

    return result_dict

def initial_SSL(dbscan_centers_df, coords, gene_expression, SSL, directions, label_dict):

    seeds = list(zip(dbscan_centers_df['center_x'], dbscan_centers_df['center_y']))

    for seed in seeds:
        distances = []
        x, y = seed

        results = calculate_initial_distances(x, y, coords, gene_expression, directions, label_dict)

        for index, neighbour_coord, distance, initial_seed in zip(results['index'], results['neighbour_coordinates'], results['distances'], results['initial_seed']):
            SSL['index'].append(index)
            SSL['neighbour_coordinates'].append(neighbour_coord)
            SSL['distances'].append(distance)
            SSL['initial_seed'].append(initial_seed)

    sorted_indices = np.argsort(SSL['distances'])

    SSL['index'] = [SSL['index'][i] for i in sorted_indices]
    SSL['distances'] = [SSL['distances'][i] for i in sorted_indices]
    SSL['neighbour_coordinates'] = [SSL['neighbour_coordinates'][i] for i in sorted_indices]
    SSL['initial_seed'] = [SSL['initial_seed'][i] for i in sorted_indices]
    
    return SSL

def get_neighbours_fast(x, y, points_array, tree, k=8):
    distances, indices = tree.query([[x, y]], k=k)
    neighbours = points_array[indices[0]]
    return neighbours

def calculate_initial_distances_kdtree(x, y, coords, gene_expression, tree, label_dict, k=8):
    result_dict = {'index': [], 'neighbour_coordinates': [], 'distances': [], 'initial_seed': []}
    initial_seed = (x, y)
    neighbours = get_neighbours_fast(x, y, coords[['x', 'y']].values, tree, k)

    for neighbour in neighbours:
        neighbour_x, neighbour_y = neighbour
        if (neighbour_x, neighbour_y) in label_dict['coordinates']:
            continue

        neighbour_index = (np.abs(coords['x'] - neighbour_x) + np.abs(coords['y'] - neighbour_y)).idxmin()
        neighbour_gene_expression = gene_expression.loc[neighbour_index]
        distance = np.linalg.norm(gene_expression.loc[coords[(coords['x'] == x) & (coords['y'] == y)].index[0]] - neighbour_gene_expression)

        result_dict['index'].append(neighbour_index)
        result_dict['neighbour_coordinates'].append((neighbour_x, neighbour_y))
        result_dict['distances'].append(distance)
        result_dict['initial_seed'].append(initial_seed)

    return result_dict

def update_SSL_kdtree(seeds, coords, gene_expression, SSL, tree, label_dict, initial_seed, k=8):
    for seed in seeds:
        x, y = seed
        results = calculate_initial_distances_kdtree(x, y, coords, gene_expression, tree, label_dict, k)

        for index, neighbour_coord, distance in zip(results['index'], results['neighbour_coordinates'], results['distances']):
            if neighbour_coord in SSL['neighbour_coordinates']:
                index_in_SSL = SSL['neighbour_coordinates'].index(neighbour_coord)
                if distance < SSL['distances'][index_in_SSL]:
                    SSL['distances'][index_in_SSL] = distance
                    SSL['index'][index_in_SSL] = index
                    SSL['initial_seed'][index_in_SSL] = initial_seed
            else:
                SSL['index'].append(index)
                SSL['neighbour_coordinates'].append(neighbour_coord)
                SSL['distances'].append(distance)
                SSL['initial_seed'].append(initial_seed)

    sorted_indices = np.argsort(SSL['distances'])
    SSL['index'] = [SSL['index'][i] for i in sorted_indices]
    SSL['distances'] = [SSL['distances'][i] for i in sorted_indices]
    SSL['neighbour_coordinates'] = [SSL['neighbour_coordinates'][i] for i in sorted_indices]
    SSL['initial_seed'] = [SSL['initial_seed'][i] for i in sorted_indices]

    return SSL

def update_label_dict(SSL, label_dict):
    
    if SSL['index'] and SSL['neighbour_coordinates']:
        if SSL['index'][0] not in label_dict['index']:
            label_dict['index'].append(SSL['index'][0])
            label_dict['coordinates'].append(SSL['neighbour_coordinates'][0])
            ssl_initial_seed = SSL['initial_seed'][0]
            index = label_dict['coordinates'].index(ssl_initial_seed)
            label_value = label_dict['label'][index]
            label_dict['label'].append(label_value)

        for key in SSL:
            if SSL[key]:
                SSL[key].pop(0)
    else:
        print("Warning: 'index' or 'neighbour_coordinates' is empty.")
    return SSL

def visualize_results(label_df, sorted_labels):
    fig, ax = plt.subplots(figsize=(8, 6))

    cmap = cm.get_cmap('tab20', len(sorted_labels))

    for i, label in enumerate(sorted_labels):
        label_data = label_df[label_df['label'] == label]
        color = cmap(i)
        ax.scatter(label_data['x'], label_data['y'], color=color, s=8, label=f'Segment {label}')
    
    ax.set_title('SRG results')
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small', markerscale=2)
    
    for spine in ax.spines.values():
        spine.set_linewidth(0)
    
    plt.tight_layout()
    plt.show()

def clustering(dbscan_centers_df, coords, st_highly_variable_genes_df):

    """
    Performs Seeded Region Growing clustering to refine initial cluster assignments.

    This function uses a seeded region growing algorithm that iteratively expands cluster labels 
    from initial DBSCAN centers. It leverages KDTree for fast neighbor searches and fills in any 
    missing labels by assigning the label of the nearest labeled coordinate.

    Parameters:
        dbscan_centers_df (pandas.DataFrame):
            DataFrame containing initial cluster centers from DBSCAN. Expected to include coordinate 
            information (e.g., as tuples in a 'coordinates' column).

        coords (pandas.DataFrame):
            DataFrame of spatial coordinates with at least 'x' and 'y' columns. Its index represents 
            unique identifiers for the spots.

        st_highly_variable_genes_df (pandas.DataFrame):
            DataFrame containing spatial transcriptomics data of highly variable genes. Its index is 
            used to determine the number of iterations for region growing.

    Returns:
        pandas.DataFrame:
            A DataFrame that contains the final clustering labels for each coordinate, with additional 
            columns for 'coordinates', 'x', 'y', and 'label'. This DataFrame is also used for visualization.
    """

    print("Running Seeded Region Growing...")
    neighbours_rules = [(-2, 0), (2, 0), (0, -2), (0, 2), (-1, -1), (1, -1), (-1, 1), (1, 1)]
    label_dict = {'coordinates': [], 'label': []}
    SSL = {'index': [], 'neighbour_coordinates': [], 'distances': [], 'initial_seed': []}

    tree = KDTree(coords[['x', 'y']].values)
    label_dict = cluster_label(dbscan_centers_df)
    SSL = initial_SSL(dbscan_centers_df, coords, st_highly_variable_genes_df, SSL, neighbours_rules, label_dict)

    for i in range(len(st_highly_variable_genes_df)):
        new_SSL = update_label_dict(SSL, label_dict)

        if new_SSL['neighbour_coordinates']:
            new_seed = [new_SSL['neighbour_coordinates'][0]]
            initial_seed = new_SSL['initial_seed'][0]
            SSL = update_SSL_kdtree(new_seed, coords, st_highly_variable_genes_df, new_SSL, tree, label_dict, initial_seed)
        else:
            print("The neighbour_coordinates list is empty. Exiting the loop.")
            break

    label_df = pd.DataFrame(label_dict, columns=['index', 'coordinates', 'label']).set_index('index')

    missing_indexes = coords.index.difference(label_df.index)
    missing_coords = coords.loc[missing_indexes, ['x', 'y']]

    min_distance_coordinates_list = []
    min_distances = []

    for index, target_point in missing_coords.iterrows():
        distances = euclidean_distances(label_df['coordinates'].tolist(), [target_point.tolist()])
        min_distance_index = distances.argmin()
        min_distance_coordinates = label_df.iloc[min_distance_index]['coordinates']

        min_distances.append(distances.min())
        min_distance_coordinates_list.append(min_distance_coordinates)

    missing_coords['min_distance_coordinates'] = min_distance_coordinates_list
    missing_coords['label'] = missing_coords['min_distance_coordinates'].apply(
        lambda coord: label_df[label_df['coordinates'] == coord]['label'].values[0]
    )

    new_coordinates = list(zip(missing_coords['x'], missing_coords['y']))
    new_labels = missing_coords['label']
    new_data = {'coordinates': new_coordinates, 'label': new_labels}
    new_df = pd.DataFrame(new_data)

    label_df = pd.concat([label_df, new_df], axis=0, ignore_index=False)
    label_df['x'], label_df['y'] = zip(*label_df['coordinates'])

    unique_labels = label_df['label'].unique()
    sorted_labels = sorted(unique_labels)

    visualize_results(label_df, sorted_labels)
    return label_df
