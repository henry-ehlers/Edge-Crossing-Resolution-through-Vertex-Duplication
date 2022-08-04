from src.graph_simulation import *
from src.graph_drawing import *
from src.vertex_splitting import *
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
    positions = embed_graph(graph=graph, embedding=embedding)
    draw_graph(graph=graph, positions=positions)

    # Save Drawing
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed)
    save_drawn_graph(output_path)

    # Check Edge Crossings
    edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
    vertex_sum = sum(vertex_crossings)
    edge_sum = 0
    for edge_crossing in edge_crossings:
        if edge_crossing is None: continue
        edge_sum += len(edge_crossing)

    print(vertex_sum) #  is four times larger than edge sum because 4 nodes map to one crossing
    print(edge_sum)