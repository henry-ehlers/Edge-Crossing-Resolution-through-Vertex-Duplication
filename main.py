from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
import networkx as nx


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
    assert vertex_edge_crossing_equality(vertex_crossings, edge_crossings), \
        "Vertex and Edge Crossing Numbers not equivalent. Sum(Cr(Vi)) / 4 = Sum(Cr(Ej))"

    # Find vertex involved in the largest number of edge crossings
    target_vertex_index = get_target_vertex_index(vertex_crossings, graph)
    print("Split Target Vertex = Vertex #{}".format(target_vertex_index))

    # Planarize Graph
    planar_graph, planar_positions = planarize_graph(graph, positions, edge_crossings)
    planar_edge_crossings, planar_vertex_crossings = locate_edge_crossings(planar_graph, planar_positions)
    # debug_edge_crossings(planar_graph, planar_edge_crossings)
    # assert is_without_edge_crossings(planar_graph, planar_positions), \
    #     "Graph is not edge-crossing free."

    # Draw and Save Non-Planar rGraph
    draw_graph(graph=graph, positions=positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=0)
    save_drawn_graph(output_path)

    # Draw and Save Planar rGraph
    draw_graph(graph=planar_graph, positions=planar_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=1)
    save_drawn_graph(output_path)
