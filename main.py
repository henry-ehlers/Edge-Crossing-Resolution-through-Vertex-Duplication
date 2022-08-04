from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
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

    # Embed Graph
    positions = embed_graph(graph=graph, embedding=embedding)

    # Check Edge Crossings
    edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
    assert vertex_edge_crossing_equality(vertex_crossings, edge_crossings)

    # Find vertex involved in the largest number of edge crossings
    target_vertex_index = get_target_vertex_index(vertex_crossings, graph)
    graph.nodes[target_vertex_index]["target"] = 1

    # Draw Graph
    draw_graph(graph=graph, positions=positions)

    # Save Drawing
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed)
    save_drawn_graph(output_path)