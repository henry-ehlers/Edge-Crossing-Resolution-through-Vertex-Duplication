import networkx as nx
import numpy as np
from pathlib import Path


def create_barabasi_albert_graph(n, m, seed=22, root_dir="data/simulated"):

    #
    file_path = f = "./{}/barabasi_albert_{}_{}_{}.csv".format(root_dir, n, m, seed)

    #
    if Path(file_path).is_file():
        data = np.loadtxt(fname=file_path, dtype=int, delimiter=",")
        graph = nx.from_numpy_matrix(A=data)
        print(graph)
    else:
        graph = nx.barabasi_albert_graph(n=n, m=m, seed=seed)
        np.savetxt(fname=file_path, X=nx.to_numpy_matrix(graph).astype(int), fmt='%i', delimiter=",")\

    #
    return graph

