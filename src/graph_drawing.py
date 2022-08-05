import matplotlib.pyplot as plt
from pathlib import Path
import networkx as nx


def embed_graph(graph, embedding, n_iter=None, seed=None):
    if embedding == "kamada_kawai":
        return nx.kamada_kawai_layout(G=graph)
    elif embedding == "fruchterman-reingold":
        return nx.spring_layout(G=graph, iterations=5000, seed=24)


def color_split_vertices(graph):
    node_color_map = []
    for vertex in graph:
        if graph.nodes[vertex]["target"] == 0:
            node_color_map.append("grey")
        else:
            node_color_map.append("red")
    return node_color_map


def color_edges(graph):
    edge_color_map = []
    for edge in graph.edges:
        edge_color_map.append("grey")
    return edge_color_map


def shape_split_vertices(graph):
    node_shape_map = []
    for vertex in graph:
        if graph.nodes[vertex]["split"] == 0:
            node_shape_map.append('o')
        else:
            node_shape_map.append('o')
    return node_shape_map


def create_output_path(embedding, n_vertices, m_edges, seed, n_splits=0):
    output_directory = "./drawings/{}/barabasi_albert_{}_{}_{}".format(embedding, n_vertices, m_edges, seed)
    Path(output_directory).mkdir(parents=True, exist_ok=True)
    file_name = "{}.png".format(n_splits)
    return "{}/{}".format(output_directory, file_name)


def save_drawn_graph(output_path):
    plt.savefig(fname=output_path, dpi=300)
    plt.clf()


def draw_graph(graph, positions):

    # Color Map for both edges and vertices
    node_color_map = color_split_vertices(graph)
    edge_color_map = color_edges(graph)
    # node_shape_map = shape_split_vertices(graph)  # CAN ONLY BE IMPLEMEMENTED USING MULTIPLE PASSES

    # Draw Graph Embedding
    nx.draw(G=graph, pos=positions, node_color=node_color_map, edge_color=edge_color_map, node_shape='o')
    nx.draw_networkx_labels(G=graph, pos=positions)
