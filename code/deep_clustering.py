from clustpy.deep import ENRC
from src.utils import kmersFasta
import numpy as np


def run(num_clusters_list, fasta_file, max_k, result_folder):
    results_json = {}
    for k in range(1, max_k + 1):
        results_json[k] = {}
        _, kmers = kmersFasta(fasta_file, k=k, transform=None, reduce=True)
        kmers_normalized = np.transpose((np.transpose(kmers) / np.linalg.norm(kmers, axis=1)))
        cluster_assignments, verbose_enrc_model, embeddings = perform_deep_clustering(kmers_normalized, num_clusters_list)
        results_json[k] = (cluster_assignments, verbose_enrc_model, embeddings)
        print(f"Finished processing k = {k}", flush=True)

    return results_json




def perform_deep_clustering(data, n_clusters_list):


    print("training started")
    # Initialize the Verbose ENRC model
    verbose_enrc_model = ENRC(n_clusters_list)

    model = verbose_enrc_model.fit(data)
    # embeddings = enrc.encoder.predict(data)

    embeddings = model.autoencoder.encode(data)
    # Get the cluster assignments with verbose output
    cluster_assignments = verbose_enrc_model.predict(data)

    embeddings = embeddings.cpu().detach().numpy()

    return cluster_assignments, verbose_enrc_model, embeddings
