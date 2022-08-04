from src.graph_simulation import *
from src.graph_drawing import *
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Define Input Parameters
    embedding = "kamada_kawai"
    n_vertices = 20
    m_edges = 3
    seed = 1

    # Create Simulated Graph
    graph = create_barabasi_albert_graph(n=n_vertices, m=m_edges, seed=seed)
    graph.nodes[0]["split"] = 1  # Arbitrarily color node as though it were split

    # Embed and Draw Graph
    draw_graph(graph=graph, embedding=embedding)

    # Save Drawing
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed)
    save_drawn_graph(output_path)

