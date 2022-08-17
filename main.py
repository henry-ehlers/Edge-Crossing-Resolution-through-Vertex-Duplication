import numpy as np
import timeit
import sys

from py2d.Math import Polygon, Vector
from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
from src.diagnostics import *
from src.line_segments import *
from src.vertex_splitting import *
from src.faces import *


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Command Line Arguments
    cmd_args = sys.argv

    # Define Input Parameters
    embedding = "kamada_kawai"
    n_vertices = int(sys.argv[1])
    m_edges = int(sys.argv[2])
    seed = int(sys.argv[3])

    # Diagnostics Files
    diagnostics_directory = "./output/diagnostics"
    diagnostics_file=f"barabasi_albert_{n_vertices}_{m_edges}_{seed}"

    P = Polygon.from_tuples([(2.00, 3.00), (3.00, 3.00), (3.00, 2.00), (4.00, 2.00), (4.00, 4.00), (2.00, 4.00)])
    P = Polygon.convex_decompose(P)
    print(P)

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
    vertex_identification_start_time = timeit.default_timer()
    target_vertex = get_target_vertex(vertex_crossings, graph)
    target_vertex_adjacency = list(graph.neighbors(target_vertex))

    # Delete Target Vertex from Graph
    remaining_edge_crossings = get_remaining_edge_crossings(graph, edge_crossings, target_vertex)
    remaining_graph, remaining_positions = remove_target_vertex(graph, positions, target_vertex)
    vertex_identification_time = timeit.default_timer() - vertex_identification_start_time
    print("Split Target Vertex = Vertex #{}".format(target_vertex))
    print("Target Adjacency    = {}".format(target_vertex_adjacency))

    # Draw and Save Graph with target removed
    initial_drawing_start_time = timeit.default_timer()
    draw_graph(graph=remaining_graph, positions=remaining_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=1)
    save_drawn_graph(output_path)
    initial_drawing_time = timeit.default_timer() - initial_drawing_start_time

    # Save time taken
    save_time(times=[vertex_identification_time, initial_drawing_time],
              labels=["identifying_target_vertex", "drawing_initial_embedding"],
              file_name=diagnostics_file,
              directory=diagnostics_directory,
              overwrite=True)

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
    planarization_start_time = timeit.default_timer()

    # Planarize Graph
    plane_graph, plane_positions, virtual_edge_set = planarize_graph(
        graph=remaining_graph, positions=remaining_positions, edge_crossings=remaining_edge_crossings)
    planar_edge_crossings, planar_vertex_crossings = locate_edge_crossings(plane_graph, plane_positions)
    planarization_time = timeit.default_timer() - planarization_start_time
    # TODO: find outer plane by finding edges that are involved in only 1 triangle -> those form the outer face bound
    # TODO: art gallery problem: split outer face into subfaces using virtual vertices
    # TODO: figure out how to deal with disconnected vertices -> do they mess up the art gallery problem?

    # Draw and Save Planar rGraph
    planar_drawing_start_time = timeit.default_timer()
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=3)
    save_drawn_graph(output_path)
    planar_drawing_time = timeit.default_timer() - planar_drawing_start_time

    # Save time taken
    save_time(times=[planarization_time, planar_drawing_time],
              labels=["planarization_of_initial_graph", "drawing_of_planar_graph"],
              file_name=diagnostics_file,
              directory=diagnostics_directory)

    # FACE IDENTIFICATION ----------------------------------------------------------------------------------------------
    face_detection_start_time = timeit.default_timer()

    # # Locate faces and best two for target face
    # TODO: list of nones in found faces
    # TODO: art gallery problem: make sure face allows for all incident vertices to be 'seen', else split up
    faces = find_all_faces(graph=plane_graph)
    face_edge_map = build_face_to_edge_map(plane_graph, faces)
    face_detection_start_time = timeit.default_timer() - face_detection_start_time
    ordered_face_edges = get_ordered_face_edges(faces, plane_graph)
    [print(f"{key} - {ordered_face_edges[key]}") for key in ordered_face_edges.keys()]
    for face in faces:
        print(f"Current {face}: {ordered_face_edges[face]}")
        decomposition = is_face_convex(face, ordered_face_edges[face], plane_positions)
        if not decomposition[0]:
            [print(decomp) for decomp in decomposition[1]]
    sys.exit()

    # Find Best Face
    face_selection_start_time = timeit.default_timer()
    face_incidences = find_face_vertex_incidence(faces, target_vertex_adjacency)
    max_incidence, selected_faces = get_maximally_incident_faces(face_incidences)
    face_selection_time = timeit.default_timer() - face_selection_start_time
    print("#Face Target Sets Found: {}".format(len(selected_faces)))
    print("Maximum Target Incidence: {}".format(max_incidence))

    # Save time taken
    save_time(times=[face_detection_start_time, face_selection_time],
              labels=["detecting_faces", "selecting_faces"],
              file_name=diagnostics_file,
              directory=diagnostics_directory)

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
    #print(f"remaining vertices: {remaining_vertex_connections}")

    # Draw and Save Planar rGraph
    face_graph_start_time = timeit.default_timer()
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=4)
    save_drawn_graph(output_path)
    face_graph_time = timeit.default_timer() - face_graph_start_time

    # Save time taken
    save_time(times=[face_graph_time],
              labels=["drawing_face_graph"],
              file_name=diagnostics_file,
              directory=diagnostics_directory)

    # All Line Segments ------------------------------------------------------------------------------------------------
    line_segment_start_time = timeit.default_timer()

    # Create line-segments between all vertices now already connected by edges or virtual edge sets
    all_segment_graph, all_segment_positions = draw_all_line_segments(
        plane_graph, plane_positions, virtual_edge_set)
    line_segment_time = timeit.default_timer() - line_segment_start_time

    line_segment_graph_start_time = timeit.default_timer()
    draw_graph(graph=all_segment_graph, positions=all_segment_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=5)
    save_drawn_graph(output_path)
    line_segment_graph_time = timeit.default_timer() - line_segment_graph_start_time

    # Save time taken
    save_time(times=[line_segment_time, line_segment_graph_time],
              labels=["projecting_all_line_segments", "drawing_line_segments"],
              file_name=diagnostics_file,
              directory=diagnostics_directory)

    # Limited Line Segments --------------------------------------------------------------------------------------------

    # Keep only segments that pass through faces of interest
    culling_start_time = timeit.default_timer()
    culled_segment_graph, culled_segment_positions, face_intersection_map = cull_all_line_segment_graph(
        all_segment_graph, all_segment_positions, selected_face_set, face_edge_map)
    culling_time = timeit.default_timer() - culling_start_time

    culled_graph_start_time = timeit.default_timer()
    draw_graph(graph=culled_segment_graph, positions=culled_segment_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=6)
    save_drawn_graph(output_path)
    culled_graph_time = timeit.default_timer() - culled_graph_start_time

    # Save time taken
    save_time(times=[culling_time, culled_graph_time],
              labels=["culling_line_segments", "drawing_culled_line_segments"],
              file_name=diagnostics_file,
              directory=diagnostics_directory)

    # Create Subfaces --------------------------------------------------------------------------------------------------

    #
    subface_creation_start_time = timeit.default_timer()
    face_graph, face_graph_positions, face_graph_virtual_edge_associations, face_vertex_map = create_subface_graph(
        culled_segment_graph, culled_segment_positions, selected_face_set, face_intersection_map)
    subface_creation_time = timeit.default_timer() - subface_creation_start_time

    draw_graph(graph=face_graph, positions=face_graph_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=7)
    save_drawn_graph(output_path)

    # Planarize the subfaces
    planarize_subfaces_start_time = timeit.default_timer()
    face_edge_crossings, face_vertex_crossings = locate_edge_crossings(face_graph, face_graph_positions)
    # TODO: list of NONE's in the plane_face_virtual_edge_set
    plane_face_graph, plane_face_graph_positions, plane_face_virtual_edge_set = planarize_graph(
        graph=face_graph, positions=face_graph_positions, edge_crossings=face_edge_crossings)
    planarize_subfaces_time = timeit.default_timer() - planarize_subfaces_start_time

    draw_graph(graph=plane_face_graph, positions=plane_face_graph_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=8)
    save_drawn_graph(output_path)

    # Save time taken
    save_time(times=[subface_creation_time, planarize_subfaces_time],
              labels=["creating_subface_graph", "planarizing_subface_graph"],
              file_name=diagnostics_file,
              directory=diagnostics_directory)

    # Determine Split Vertex Placement ---------------------------------------------------------------------------------

    # Get all subfaces of all target faces
    face_subface_start_time = timeit.default_timer()
    plane_graph_sub_faces = find_all_subfaces(plane_face_graph, plane_face_virtual_edge_set, face_vertex_map)
    face_subface_time = timeit.default_timer() - face_subface_start_time

    # Obtain all remaining vertices to be connected
    subface_centroid_start_time = timeit.default_timer()
    subface_centroids = get_split_vertex_locations(plane_face_graph_positions, plane_graph_sub_faces)
    subface_centroid_time = timeit.default_timer() - subface_centroid_start_time

    # Calculate induced edge crossings using the graph with the deleted vertex
    induced_crossings_start_time = timeit.default_timer()
    induced_edge_crossings = calculate_induced_edge_crossings(
        remaining_graph, remaining_positions, subface_centroids, remaining_vertex_connections)
    pair_induced_crossings, neighbor_assignment = get_split_vertex_pairs(induced_edge_crossings)
    induced_crossings_time = timeit.default_timer() - induced_crossings_start_time

    # Find vertex pair which minimizes joint induced edge crossings
    vertex_selection_start_time = timeit.default_timer()
    min_induced_crossings, target_subfaces, ties = select_vertex_splits(pair_induced_crossings)
    vertex_selection_time = timeit.default_timer() - vertex_selection_start_time
    print(f"Minimum Induced Crossings: {min_induced_crossings}")
    print(f"target Subfaces: {target_subfaces}")
    print(f"Number of ties: {ties}")

    # Save time taken
    save_time(times=[face_subface_time, subface_centroid_time, induced_crossings_time, vertex_selection_time],
              labels=["linking_faces_to_subfaces", "calculating_subface_centroids",
                      "calculcating_induced_edge_crossings", "selecting_target_subfaces"],
              file_name=diagnostics_file,
              directory=diagnostics_directory)

    # Place Split Vertices ---------------------------------------------------------------------------------------------

    # TODO: add function "melt_frozen_set" to add stuff to it

    split_graph, splot_positions, connected_vertices = place_split_vertices(
        remaining_graph, remaining_positions, target_vertex, target_vertex_adjacency, remaining_vertex_connections,
        neighbor_assignment, subface_centroids, target_subfaces, plane_graph_sub_faces)

    draw_graph(graph=split_graph, positions=splot_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=9)
    save_drawn_graph(output_path)
