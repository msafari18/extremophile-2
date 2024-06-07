import torch
from torchvision import datasets, transforms
import numpy as np
from clustpy.deep import ENRC
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

# Load the MNIST dataset
transform = transforms.Compose([transforms.ToTensor(), transforms.Normalize((0.5,), (0.5,))])
mnist_data = datasets.MNIST(root='./data', train=True, download=True, transform=transform)

# Convert the dataset to a NumPy array
data = mnist_data.data.numpy().reshape(-1, 28*28) / 255.0  # Normalize pixel values to [0, 1]
data = torch.tensor(data, dtype=torch.float)

# Print the shape of the data
print("Dataset shape:", data.shape)

# Define the number of clusters for each aspect
n_clusters_list = [2, 10, 2]  # Clusters for digit shape/style, digit class, and intensity/stroke width

# Initialize the ENRC model
enrc_model = ENRC(n_clusters_list)

# Fit the model
enrc_model.fit(data)

# Get the cluster assignments
cluster_assignments = enrc_model.predict(data)

# Get the embeddings
embeddings = enrc_model.embeddings_

# Print the shape of the embeddings to verify
print("Embeddings shape:", embeddings.shape)

# Apply t-SNE to the embeddings
tsne = TSNE(n_components=2)
embeddings_2d = tsne.fit_transform(embeddings)

# Plot the embeddings (using the first aspect for coloring)
plt.figure(figsize=(10, 7))
plt.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=cluster_assignments[0], cmap='viridis')
plt.colorbar(label='Cluster Label')
plt.title("t-SNE visualization of embeddings")
plt.xlabel("t-SNE Component 1")
plt.ylabel("t-SNE Component 2")
plt.show()
