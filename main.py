from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
import networkx as nx
import numpy as np


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

    # Collect and Check Edge Crossings
    edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
    debug_edge_crossings(graph, edge_crossings, positions)
    assert vertex_edge_crossing_equality(vertex_crossings, edge_crossings), \
        "Vertex and Edge Crossing Numbers not equivalent. Sum(Cr(Vi)) / 4 = Sum(Cr(Ej))"

    # Check rectangular bounds of drawing
    min_coord, max_coord = find_embedding_rectangle(graph, positions)
    print("Minimum Coordinates: {}".format(min_coord))
    print("Maximum Coordinates: {}".format(max_coord))

    # Draw and Save Non-Planar rGraph
    draw_graph(graph=graph, positions=positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=0)
    save_drawn_graph(output_path)

    # Find vertex involved in the largest number of edge crossings
    target_vertex_index = get_target_vertex_index(vertex_crossings, graph)
    print("Split Target Vertex = Vertex #{}".format(target_vertex_index))

    # Planarize Graph
    plane_graph, plane_positions = planarize_graph(graph, positions, edge_crossings)
    planar_edge_crossings, planar_vertex_crossings = locate_edge_crossings(plane_graph, plane_positions)
    debug_edge_crossings(plane_graph, planar_edge_crossings, plane_positions)
    is_planar, plane_graph = nx.check_planarity(plane_graph)
    find_faces(plane_graph, plane_positions)

    # Check rectangular bounds of drawing
    min_coord, max_coord = find_embedding_rectangle(plane_graph, plane_positions)
    print("Minimum Coordinates: {}".format(min_coord))
    print("Maximum Coordinates: {}".format(max_coord))

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=1)
    save_drawn_graph(output_path)
