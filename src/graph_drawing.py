import matplotlib.pyplot as plt
from pathlib import Path
import networkx as nx


def get_embedding_square(graph, positions, scaler=1):
    maximum = float("-inf")
    for vertex in graph.nodes:
        largest_node_coordinate = max([abs(x) for x in positions[vertex]])
        if largest_node_coordinate <= maximum:
            continue
        maximum = largest_node_coordinate
    maximum *= scaler
    return (-maximum, -maximum), (-maximum, maximum), (maximum, maximum), (maximum, -maximum)


def find_embedding_rectangle(graph, positions):

    # Initialize minimum and maximum coordinates using first (arbitrary) selected vertex
    min_coordinates, max_coordinates = [positions[0][0], positions[0][1]], [positions[0][0], positions[0][1]]

    # Iterate over all vertices to find maximum and minimum x/y values
    for vertex_index in graph.nodes:

        # Check minimum and maximum x-coordinate values
        if positions[vertex_index][0] < min_coordinates[0]:
            min_coordinates[0] = positions[vertex_index][0]
        elif positions[vertex_index][0] > max_coordinates[0]:
            max_coordinates[0] = positions[vertex_index][0]

        # Check minimum and maximum y-coordinate values
        if positions[vertex_index][1] < min_coordinates[1]:
            min_coordinates[1] = positions[vertex_index][1]
        elif positions[vertex_index][1] > max_coordinates[1]:
            max_coordinates[1] = positions[vertex_index][1]

    # Return maximum and minimum
    return min_coordinates, max_coordinates


def embed_graph(graph, embedding, n_iter=None, seed=None):
    if embedding == "kamada_kawai":
        return nx.kamada_kawai_layout(G=graph)
    elif embedding == "fruchterman-reingold":
        return nx.spring_layout(G=graph, iterations=n_iter, seed=seed)


def get_graph_entity_data(entity_dictionary, entity_key, data_key, default_value):

    # Ensure the entity exists in dictionary
    assert entity_dictionary.get(entity_key, None) is not None, \
        f"Entity {entity_key} does not exist in provided dictionary."

    # Try to access the entity's request data
    try:
        data = entity_dictionary[entity_key][data_key]
    except KeyError:
        data = default_value

    # Return Retrieved or Default Values
    return data


def color_split_vertices(graph):
    node_color_map = []
    for vertex in graph:
        if get_graph_entity_data(graph.nodes, vertex, "split", 0) == 1:
            color = "red"
        elif get_graph_entity_data(graph.nodes, vertex, "boundary", 0) == 1:
            color = "blue"
        elif get_graph_entity_data(graph.nodes, vertex, "virtual", 0) == 1:
            color = "lightgrey"
        else:
            color = "lightgrey"
        node_color_map.append(color)
    return node_color_map


def color_edges(graph):
    edge_color_map = []
    for edge in graph.edges:
        if get_graph_entity_data(graph.edges, edge, "segment", 0) == 1:
            color = "green"
        elif get_graph_entity_data(graph.edges, edge, "target", 0) == 1:
            color = "blue"
        elif get_graph_entity_data(graph.edges, edge, "virtual", 0) == 0:
            color = "black"
        elif get_graph_entity_data(graph.edges, edge, "virtual", 0) == 1:
            color = "lightgrey"
        else:
            color = "lightgrey"
        edge_color_map.append(color)
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


def draw_graph(graph, positions, labels=None):

    # Color Map for both edges and vertices
    node_color_map = color_split_vertices(graph)
    edge_color_map = color_edges(graph)
    # node_shape_map = shape_split_vertices(graph)  # CAN ONLY BE IMPLEMEMENTED USING MULTIPLE PASSES

    # Draw Graph Embedding
    plt.figure(3, figsize=(30, 30))
    nx.draw(G=graph, pos=positions, node_color=node_color_map, edge_color=edge_color_map, node_shape='o', node_size=75)
    if labels is not None:
        nx.draw_networkx_labels(G=graph, pos=positions, labels=labels, font_size=15)
    else:
        nx.draw_networkx_labels(G=graph, pos=positions, font_size=15)
