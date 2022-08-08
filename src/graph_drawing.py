import matplotlib.pyplot as plt
from pathlib import Path
import networkx as nx


def find_faces(plane_graph, plane_positions):

    # Ensure the provided graph is planar
    assert nx.is_planar(plane_graph), \
        "Provided Graph is not planar."

    # Initialize set of all faces (sets of vertices)
    faces = set()

    # Iterate over each vertex in the drawing
    for origin_vertex in range(0, plane_graph.number_of_nodes()):

        # Collect the neighbors [w] of a vertex [v] in clock-wise order
        cw_neighbors = list(plane_graph.neighbors_cw_order(origin_vertex))

        # For each clockwise-sorted neighbor w, find the face to the right of the half-edge (v,w)
        for adjacent_vertex in cw_neighbors:
            face = frozenset(plane_graph.traverse_face(origin_vertex, adjacent_vertex))
            faces.add(face)

    # Return set of faces (frozen sets of vertices)
    return faces


def find_embedding_rectangle(graph, positions):

    # Initialize minimum and maximum coordinates using first (arbitrary) selected vertex
    min_coordinates, max_coordinates = [positions[0][0], positions[0][1]], [positions[0][0], positions[0][1]]
    print("Number of Nodes: {}".format(len(graph.nodes)))

    # Iterate over all vertices to find maximum and minimum x/y values
    for vertex_index in range(1, len(graph.nodes)):

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
        return nx.spring_layout(G=graph, iterations=5000, seed=24)


def color_split_vertices(graph):
    node_color_map = []
    for vertex in graph:
        if graph.nodes[vertex]["virtual"] == 0 and graph.nodes[vertex]["target"] == 0:
            color = "black"
        elif graph.nodes[vertex]["target"] == 1:
            color = "red"
        else:
            color = "lightgrey"
        node_color_map.append(color)
    return node_color_map


def color_edges(graph):
    virtual = nx.get_edge_attributes(graph, "virtual")
    edge_color_map = []
    for edge in graph.edges:
        if virtual[edge] == 0:
            color = "black"
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


def draw_graph(graph, positions):

    # Color Map for both edges and vertices
    node_color_map = color_split_vertices(graph)
    edge_color_map = color_edges(graph)
    # node_shape_map = shape_split_vertices(graph)  # CAN ONLY BE IMPLEMEMENTED USING MULTIPLE PASSES

    # Draw Graph Embedding
    plt.figure(3, figsize=(20, 20))
    nx.draw(G=graph, pos=positions, node_color=node_color_map, edge_color=edge_color_map, node_shape='o', node_size=7)
    nx.draw_networkx_labels(G=graph, pos=positions, font_size=10)
