import matplotlib.pyplot as plt
import networkx as nx
import random as rnd
import numpy as np


def create_barabasi_albert_graph(n, m, seed=22, root_dir="data/simulated"):
    file_path = f = "./{}/barabasi_albert_{}_{}_{}.csv".format(root_dir, n, m, seed)
    graph = nx.barabasi_albert_graph(n=n, m=m, seed=seed)
    np.savetxt(fname=file_path, X=nx.to_numpy_matrix(graph).astype(int), fmt='%i', delimiter=",")
    return graph

