from src.create_random_graphs import create_barabasi_albert_graph
import matplotlib.pyplot as plt
import networkx as nx

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    graph = create_barabasi_albert_graph(n=20, m=3, seed=1)
    print("here")
    #pos = nx.spring_layout(G=graph, iterations=5000, seed=24)
    pos = nx.kamada_kawai_layout(G=graph)

    # Color Map for both edges and vertices
    node_color_map = []
    edge_color_map = []
    for vertex in graph:
        node_color_map.append("grey")
    for edge in graph.edges:
        edge_color_map.append("grey")

    # Plot Graph
    nx.draw(G=graph, pos=pos, node_color=node_color_map, edge_color=edge_color_map)
    nx.draw_networkx_labels(G=graph, pos=pos)
    plt.savefig("pre_crossing.png")
    plt.clf()

