import matplotlib.pyplot as plt
import networkx as nx
import random as rnd
import numpy as np


def create_barabasi_albert_graph(n, m, root_dir="data/simulated"):
    file_path = f = "./{}/barabasi_albert_{}_{}.csv".format(root_dir, n, m)
    graph = nx.barabasi_albert_graph(n=n, m=m, seed=22)
    np.savetxt(fname=f, X=nx.to_numpy_matrix(G).astype(int), fmt='%i', delimiter=",")
