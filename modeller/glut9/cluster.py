#!/usr/bin/env python

# Script to calculate RMSD between multiple structural models
# in a directory using BioPython's Superimposer class

import os
from Bio import PDB
from Bio.PDB import Superimposer
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from collections import Counter
import matplotlib.pyplot as plt
import shutil
# Assuming rmsd_matrix is your RMSD matrix

def calculate_rmsd(model1, model2):
    # Create a Superimposer object
    super_imposer = Superimposer()

    # Get the atom structures of the models
    atoms1 = [atom for atom in model1.get_atoms()if atom.id == 'CA']
    atoms2 = [atom for atom in model2.get_atoms()if atom.id == 'CA']

    # Set the atom coordinates for superimposition
    super_imposer.set_atoms(atoms1, atoms2)

    # Calculate RMSD
    rmsd = super_imposer.rms

    return rmsd

def calculate_rmsd_matrix(models):
    num_models = len(models)
    rmsd_matrix = [[0.0] * num_models for _ in range(num_models)]

    for i in range(num_models):
        for j in range(i + 1, num_models):
            rmsd = calculate_rmsd(models[i], models[j])
            rmsd_matrix[i][j] = rmsd
            rmsd_matrix[j][i] = rmsd  # The matrix is symmetric

    return rmsd_matrix

# Specify the path to your directory containing the structural models
directory_path = '/home/sree/git/integrative-modeling-trials/modeller/glut9'
output_directory = '/home/sree/git/integrative-modeling-trials/modeller/glut9/output-structures'

# List all files in the directory that start with 'glut9.B'
file_list = [f for f in os.listdir(directory_path) if f.startswith('glut9.B') and os.path.isfile(os.path.join(directory_path, f))]

# Create full paths for each filele
file_paths = [os.path.join(directory_path, file) for file in file_list]
# Parse structures from files
structures = []
pdb_parser = PDB.PDBParser(QUIET=True)
for file_path in file_paths:
    structure = pdb_parser.get_structure(os.path.splitext(os.path.basename(file_path))[0], file_path)
    structures.append(structure)

# Calculate RMSD matrix
rmsd_matrix = calculate_rmsd_matrix(structures)

# Print or further analyze the rmsd_matrix as needed
print(rmsd_matrix)

# Convert the RMSD matrix to a condensed distance matrix
condensed_distance = np.array([rmsd_matrix[i][j] for i in range(len(rmsd_matrix)) for j in range(i + 1, len(rmsd_matrix[i]))])

# Perform hierarchical clustering
linkage_matrix = linkage(condensed_distance, method='ward')

# Set a threshold for clustering
threshold = 2

# Cut the dendrogram to obtain clusters
clusters = fcluster(linkage_matrix, threshold, criterion='distance')

# Count the occurrences of each cluster
cluster_counts = Counter(clusters)

# Visualize the dendrogram
plt.figure(figsize=(12, 6))
dendrogram(linkage_matrix, color_threshold=threshold, above_threshold_color='gray')
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Model Index')
plt.ylabel('Distance')
plt.savefig('inward_closed.png', dpi =300)
plt.show()

# Print the cluster assignments and counts
print("Cluster Assignments:", clusters)
print("Unique Clusters and Counts:", cluster_counts)

input_directory = '/home/sree/git/integrative-modeling-trials/modeller/glut9'
output_directory = '/home/sree/git/integrative-modeling-trials/modeller/glut9/representative_models'

# Ensure the output directory exists, create it if necessary
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# List all files in the directory that start with 'glut9.B'
file_list = [f for f in os.listdir(directory_path) if f.startswith('glut9.B') and os.path.isfile(os.path.join(directory_path, f))]

# Calculate the centroid or medoid of each cluster
centroids = []
for cluster_label in set(clusters):
    cluster_indices = np.where(clusters == cluster_label)[0]
    cluster_rmsds = [rmsd_matrix[i][j] for i in cluster_indices for j in cluster_indices if i < j]
    centroid_index = cluster_indices[np.argmin(np.mean(cluster_rmsds, axis=0))]
    centroids.append(centroid_index)

# Choose models based on a threshold RMSD
threshold_rmsd = 6.0
representative_models = [(file_paths[centroid], rmsd_matrix[centroid]) for centroid in centroids if np.min(rmsd_matrix[centroid]) < threshold_rmsd]

# Save representative models to the output directory
for file_path, _ in representative_models:
    filename = os.path.basename(file_path)
    output_path = os.path.join(output_directory, filename)
    shutil.copy(file_path, output_path)

# Print the names and RMSD values of representative models
for file_path, rmsd_values in representative_models:
    print(f"Representative Model: {file_path}, RMSD: {rmsd_values}")