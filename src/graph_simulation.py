import networkx as nx
import numpy as np
from pathlib import Path


def create_barabasi_albert_graph(n, m, seed, root_dir="data/simulated", type="barabasi_albert"):

    # Define the input/output file path of the specified graph
    # file_path = f = "./{}/barabasi_albert_{}_{}_{}.csv".format(root_dir, n, m, seed)
    file_path = f = "./{}/{}_{}_{}_{}.csv".format(root_dir, type, n, m, seed)

    # If the path already exists, load it from the existing file
    if Path(file_path).is_file():

        # Load data and convert matrix into graph object
        data = np.loadtxt(fname=file_path, dtype=int, delimiter=",")
        graph = nx.from_numpy_matrix(A=data)

    # If the file does not exist, then create the graph object first, then save it to file
    else:

        # Create graph and save to file
        if type == "barabasi_albert":
                graph = nx.barabasi_albert_graph(n=n, m=m, seed=seed)
        elif type == "watts_strogatz":
            graph = nx.watts_strogatz_graph(n=n, k=m, p=0.8, seed=seed)
        np.savetxt(fname=file_path, X=nx.to_numpy_matrix(graph).astype(int), fmt='%i', delimiter=",")

    # Set Edge Attributes
    nx.set_edge_attributes(graph, 1, "real")
    nx.set_node_attributes(graph, 1, "real")

    # Return graph object
    return graph

