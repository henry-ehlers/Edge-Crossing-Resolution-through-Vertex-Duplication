import networkx as nx
import numpy as np
from pathlib import Path


def create_barabasi_albert_graph(n, m, seed, root_dir="data/simulated"):

    # Define the input/output file path of the specified graph
    file_path = f = "./{}/barabasi_albert_{}_{}_{}.csv".format(root_dir, n, m, seed)

    # If the path already exists, load it from the existing file
    if Path(file_path).is_file():

        # Load data and convert matrix into graph object
        data = np.loadtxt(fname=file_path, dtype=int, delimiter=",")
        graph = nx.from_numpy_matrix(A=data)

    # If the file does not exist, then create the graph object first, then save it to file
    else:

        # Create graph and save to file
        graph = nx.barabasi_albert_graph(n=n, m=m, seed=seed)
        np.savetxt(fname=file_path, X=nx.to_numpy_matrix(graph).astype(int), fmt='%i', delimiter=",")

    # Set Edge Attributes
    nx.set_edge_attributes(graph, 0, "virtual")
    nx.set_edge_attributes(graph, 0, "target")

    # Append the 'split' field to the graph object
    for vertex in graph:
        graph.nodes[vertex]["split"] = 0
        graph.nodes[vertex]["target"] = 0
        graph.nodes[vertex]["virtual"] = 0

    # Return graph object
    return graph

