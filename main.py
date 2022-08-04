from src.create_random_graphs import create_barabasi_albert_graph
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path


def calculate_vertex_positions(graph, embedding, n_iter=None, seed=None):
    if embedding == "kamada_kawai":
        return nx.kamada_kawai_layout(G=graph)
    elif embedding == "fruchterman-reingold":
        return nx.spring_layout(G=graph, iterations=5000, seed=24)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    embedding = "kamada_kawai"
    n_vertices = 20
    m_edges = 3
    seed = 1

    # Create Simulated Graph
    graph = create_barabasi_albert_graph(n=n_vertices, m=m_edges, seed=seed)

    # Embed Graph in 2D space to obtain vertex positions
    pos = calculate_vertex_positions(graph=graph, embedding=embedding)

    # Color Map for both edges and vertices
    node_color_map = []
    edge_color_map = []
    for vertex in graph:
        node_color_map.append("grey")
    for edge in graph.edges:
        edge_color_map.append("grey")

    # Name Currently Embedded/Drawn Graph
    n_split_events = 0
    output_directory = "./drawings/{}/barabasi_albert_{}_{}_{}".format(embedding, n_vertices, m_edges, seed)
    Path(output_directory).mkdir(parents=True, exist_ok=True)
    file_name = "{}.png".format(n_split_events)
    output_path = "{}/{}".format(output_directory, file_name)

    # Draw Graph Embedding
    nx.draw(G=graph, pos=pos, node_color=node_color_map, edge_color=edge_color_map)
    nx.draw_networkx_labels(G=graph, pos=pos)

    # Save Drawing
    plt.savefig(fname=output_path, dpi=300)
    plt.clf()

