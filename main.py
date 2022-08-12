from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
from src.line_segments import *
from src.vertex_splitting import *
from src.faces import *

import networkx as nx
import numpy as np
import sys

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Define Input Parameters
    embedding = "kamada_kawai"
    n_vertices = 10
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
    plane_graph, plane_positions, virtual_edge_set = planarize_graph(
        graph=remaining_graph, positions=remaining_positions, edge_crossings=remaining_edge_crossings)
    planar_edge_crossings, planar_vertex_crossings = locate_edge_crossings(plane_graph, plane_positions)

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=3)
    save_drawn_graph(output_path)

    # FACE IDENTIFICATION ----------------------------------------------------------------------------------------------

    # # Locate faces and best two for target face
    # TODO: list of nones in found faces
    faces = find_all_faces(graph=plane_graph)
    face_edge_map = build_face_to_edge_map(plane_graph, faces)
    [print(face) for face in faces]

    # Find Best Face
    face_incidences = find_face_vertex_incidence(faces, target_vertex_adjacency)
    max_incidence, selected_faces = get_maximally_incident_faces(face_incidences)
    print("#Face Target Sets Found: {}".format(len(selected_faces)))
    print("Maximum Target Incidence: {}".format(max_incidence))

    # FACE SELECTION ---------------------------------------------------------------------------------------------------

    # TODO: either select a single such face-set or try all that are tied
    selected_face_set_index = 0
    selected_face_set = selected_faces[selected_face_set_index]
    color_selected_faces(plane_graph, selected_face_set, face_edge_map)
    found_targets = list((selected_face_set[0].union(selected_face_set[1])).intersection(set(target_vertex_adjacency)))
    print("Selected Faces: {}".format(selected_faces[selected_face_set_index]))
    print("Incident Target Vertices: {}".format(found_targets))

    # How many neighbors have we already reached
    remaining_vertex_connections = set(target_vertex_adjacency) - set(found_targets)
    print(f"remaining vertices: {remaining_vertex_connections}")

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=4)
    save_drawn_graph(output_path)

    # All Line Segments ------------------------------------------------------------------------------------------------

    # Create line-segments between all vertices now already connected by edges or virtual edge sets
    all_segment_graph, all_segment_positions = draw_all_line_segments(
        plane_graph, plane_positions, virtual_edge_set)

    draw_graph(graph=all_segment_graph, positions=all_segment_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=5)
    save_drawn_graph(output_path)

    # Limited Line Segments --------------------------------------------------------------------------------------------

    # Keep only segments that pass through faces of interest
    culled_segment_graph, culled_segment_positions, face_intersection_map = cull_all_line_segment_graph(
        all_segment_graph, all_segment_positions, selected_face_set, face_edge_map)

    draw_graph(graph=culled_segment_graph, positions=culled_segment_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=6)
    save_drawn_graph(output_path)

    # Create Subfaces --------------------------------------------------------------------------------------------------

    #
    face_graph, face_graph_positions, face_graph_virtual_edge_associations, face_vertex_map = create_subface_graph(
        culled_segment_graph, culled_segment_positions, selected_face_set, face_intersection_map)

    draw_graph(graph=face_graph, positions=face_graph_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=7)
    save_drawn_graph(output_path)

    # Planarize the subfaces
    face_edge_crossings, face_vertex_crossings = locate_edge_crossings(face_graph, face_graph_positions)
    # TODO: list of NONE's in the plane_face_virtual_edge_set
    plane_face_graph, plane_face_graph_positions, plane_face_virtual_edge_set = planarize_graph(
        graph=face_graph, positions=face_graph_positions, edge_crossings=face_edge_crossings)

    draw_graph(graph=plane_face_graph, positions=plane_face_graph_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=8)
    save_drawn_graph(output_path)

    # Split Vertex Placement -------------------------------------------------------------------------------------------

    # Get all subfaces of all target faces
    plane_graph_sub_faces = find_all_subfaces(plane_face_graph, plane_face_virtual_edge_set, face_vertex_map)
    [print(f"{key} - {plane_graph_sub_faces[key]}") for key in plane_graph_sub_faces.keys()]

    # Obtain all remaining vertices to be connected
    subface_centroids = get_split_vertex_locations(plane_face_graph_positions, plane_graph_sub_faces)

    # Calculate induced edge crossings using the graph with the deleted vertex
    induced_edge_crossings = calculate_induced_edge_crossings(
        remaining_graph, remaining_positions, subface_centroids, remaining_vertex_connections)
    pair_induced_crossings, neighbor_assignment = get_split_vertex_pairs(induced_edge_crossings)

    # Find best vertex pair