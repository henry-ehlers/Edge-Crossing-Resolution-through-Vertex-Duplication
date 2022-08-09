from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
from src.line_segments import *
from src.faces import *

import networkx as nx
import numpy as np

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Define Input Parameters
    embedding = "kamada_kawai"
    n_vertices = 20
    m_edges = 3
    seed = 1

    # EMBED INPUT GRAPH ------------------------------------------------------------------------------------------------

    # Create Simulated Graph
    graph = create_barabasi_albert_graph(n=n_vertices, m=m_edges, seed=seed)
    positions = embed_graph(graph=graph, embedding=embedding)
    virtual_index_start = graph.number_of_nodes()

    # Collect and Check Edge Crossings
    edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
    assert vertex_edge_crossing_equality(vertex_crossings, edge_crossings), \
        "Vertex and Edge Crossing Numbers not equivalent. Sum(Cr(Vi)) / 4 = Sum(Cr(Ej))"

    # Draw and Save Non-Planar rGraph
    draw_graph(graph=graph, positions=positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=0)
    save_drawn_graph(output_path)

    # LOCATE AND REMOVE TARGET VERTEX ----------------------------------------------------------------------------------

    # Find vertex involved in the largest number of edge crossings
    target_vertex = get_target_vertex(vertex_crossings, graph)
    target_vertex_adjacency = list(graph.neighbors(target_vertex))
    print("Split Target Vertex = Vertex #{}".format(target_vertex))
    print("Target Adjacency    = {}".format(target_vertex_adjacency))

    # Delete Target Vertex from Graph
    remaining_edge_crossings = get_remaining_edge_crossings(graph, edge_crossings, target_vertex)
    remaining_graph, remaining_positions = remove_target_vertex(graph, positions, target_vertex)

    # Draw and Save Graph with target removed
    draw_graph(graph=remaining_graph, positions=remaining_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=1)
    save_drawn_graph(output_path)

    # DEBUG ------------------------------------------------------------------------------------------------------------

    # Add Nodes for each edge crossing
    debug_graph = copy.deepcopy(remaining_graph)
    debug_positions = copy.deepcopy(remaining_positions)
    n = virtual_index_start
    for edge_a in remaining_edge_crossings.keys():
        for edge_b in remaining_edge_crossings[edge_a].keys():
            debug_graph.add_node(node_for_adding=n, split=0, target=1, virtual=0)
            debug_positions[n] = np.asarray(remaining_edge_crossings[edge_a][edge_b])
            n += 1

    # Draw and Save Graph with target removed
    draw_graph(graph=debug_graph, positions=debug_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=2)
    save_drawn_graph(output_path)

    # PLANARIZATION ----------------------------------------------------------------------------------------------------

    # Planarize Graph
    plane_graph, plane_positions = planarize_graph(graph=remaining_graph,
                                                   positions=remaining_positions,
                                                   edge_crossings=remaining_edge_crossings,
                                                   starting_index=virtual_index_start)
    planar_edge_crossings, planar_vertex_crossings = locate_edge_crossings(plane_graph, plane_positions)

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=3)
    save_drawn_graph(output_path)

    # FACE IDENTIFICATION ----------------------------------------------------------------------------------------------

    # # Locate faces and best two for target face
    faces = find_all_faces(graph=plane_graph)
    face_edge_map = build_face_to_edge_map(plane_graph, faces)

    # Find Best Face
    face_incidences = find_face_vertex_incidence(faces, target_vertex_adjacency)
    max_incidence, selected_faces = get_maximally_incident_faces(face_incidences)
    print("#Face Target Sets Found: {}".format(len(selected_faces)))
    print("Maximum Target Incidence: {}".format(max_incidence))

    # Arbitrary Selection of target face-set
    # TODO: either select a single such face-set or try all that are tied
    selected_face_set_index = 0
    selected_face_set = selected_faces[selected_face_set_index]
    color_selected_faces(plane_graph, selected_face_set, face_edge_map)
    found_targets = list((selected_face_set[0].union(selected_face_set[1])).intersection(set(target_vertex_adjacency)))
    print("Selected Faces: {}".format(selected_faces[selected_face_set_index]))
    print("Incident Target Vertices: {}".format(found_targets))

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=4)
    save_drawn_graph(output_path)

    for vertex in list(plane_graph.nodes()):
        print("Vertex {}: {}".format(vertex, plane_positions[vertex]))
    # All Line Segments ------------------------------------------------------------------------------------------------
    # TODO: only connect REAL (NOT VIRTUAL) and non-edge-connected nodes to avoid overlapping lines/edges?
    # TODO: figure out how to best know which point connects to which new points below
    drawing_bounds = [-1, -1, 1, 1]
    position_1 = plane_positions[16]
    position_2 = plane_positions[17]
    print("points {} and {}".format(position_1, position_2))
    test = extend_line_to_bounds(position_1, position_2)
    print(test)
    test2 = extend_line(-1, -1, 1, 1, position_1[0], position_1[1], position_2[0], position_2[1])
    print(test2)
    start_index = plane_graph.number_of_nodes() + 1
    for index, vertex_index in enumerate(range(start_index, start_index+2)):
        print("Index {} for New Vertex {}".format(index, vertex_index))
        plane_graph.add_node(node_for_adding=vertex_index, split=0, target=0, virtual=0, boundary=1)
        if index == 0:
            plane_positions[vertex_index] = np.asarray((test[0], test[1]))
        elif index == 1:
            plane_positions[vertex_index] = np.asarray((test[2], test[3]))

    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=5)
    save_drawn_graph(output_path)

