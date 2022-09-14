import matplotlib.pyplot as plt
import networkx as nx
import itertools as it
import numpy as np
import timeit
import sys

from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
from src.diagnostics import *
from src.line_segments import *
from src.vertex_splitting import *
from src.sight_cells import *
from src.faces import *


def get_sight_cells(face_selection, target_vertices,
                    inner_faces, outer_faces,
                    sorted_inner_face_edges, sorted_outer_face_edges,
                    inner_cycles, outer_cycles,
                    graph, positions):

    # Initialize the objects in which to store each face's selected cells, their graphs, and positions
    selected_faces = [selected_face[1] for selected_face in face_selection]
    cell_selections = {face: {} for face in selected_faces}
    selected_cell_graph = {face: None for face in selected_faces}
    selected_cell_positions = {face: None for face in selected_faces}

    # Check for convexity
    for selected_face in face_selection:
        is_outer, face = selected_face[0], selected_face[1]

        # Get Sight Cells of Inner Face
        if not is_outer:
            face_edges = sorted_inner_face_edges[face]
            if is_inner_face_convex(face_edges, p_positions):
                cell_selection, cell_graph, cell_positions = {face: face.intersection(target_vertices)}, None, None
            else:
                cell_selection, cell_graph, cell_positions = None, None, None

        # Get Sight Cells of Outer Face
        else:
            cell_selection, cell_graph, cell_positions = select_outer_face_sight_cells(
                outer_faces=outer_faces,
                sorted_outer_edges=sorted_outer_face_edges,
                is_cycle=outer_cycles,
                target_vertices=target_vertices,
                graph=graph,
                positions=positions)

        # Store Selected Results
        selected_cell_graph[face] = cell_graph
        selected_cell_positions[face] = cell_positions
        cell_selections[face].update(cell_selection)

    print(f"selected_cell_graphs: {selected_cell_graph}")
    print(f"selected_cell_positions: {selected_cell_positions}")
    print(f"cell_selections: {cell_selections}")

    return cell_selections, selected_cell_graph, selected_cell_positions


def select_embedding_faces(p_graph, p_positions, target_vertices):

    # Identify the graph's faces and their incidences
    inner_faces = find_inner_faces(graph=p_graph, positions=p_positions)
    inner_face_edge_map = build_face_to_edge_map(p_graph, inner_faces)
    inner_face_incidences = find_face_vertex_incidence(inner_faces, target_vertices)
    sorted_inner_face_edges = get_ordered_face_edges(inner_faces, p_graph)

    # Find Outer Face
    outer_faces, outer_face_edges, is_cycle = find_outer_face(sorted_inner_face_edges, p_graph, p_positions)

    #
    outer_face_identifier = frozenset(set.union(*[set(outer_face) for outer_face in outer_faces]))
    sorted_outer_face__edges = {outer_face: sort_face_edges(outer_face_edges[outer_face])
                                for outer_face in outer_faces if is_cycle[outer_face]}

    # Get the incidences of the outer face
    outer_face_incidences = find_outer_face_vertex_incidence(outer_face_identifier, inner_faces, target_vertices)
    # TODO: organize data structure for outer and inner incidences
    # TODO: rank said data structure
    # TODO: figure out who to pass on the sight cell part -> needs to know in what direction to extend sight lines

    # Combine face incidences and select top set
    face_incidence_table = get_face_incidence_table(inner_face_incidences=inner_face_incidences,
                                                    outer_face_incidences=outer_face_incidences)
    print(face_incidence_table)

    # selected_faces = get_maximally_incident_faces(face_incidence_table)
    artificial_selection = face_incidence_table.loc[1]
    selected_faces = (artificial_selection.at["outer"], artificial_selection.at["face_a"]),\
                     (False, artificial_selection.at["face_b"])
    print(f"selected faces: {selected_faces}")

    cell_selections, selected_cell_graph, selected_cell_positions = get_sight_cells(
        face_selection=selected_faces,
        target_vertices=target_vertices,
        inner_faces=inner_faces,
        outer_faces=outer_faces,
        sorted_inner_face_edges=sorted_inner_face_edges,
        sorted_outer_face_edges=sorted_outer_face__edges,
        inner_cycles=None,
        outer_cycles=is_cycle,
        graph=p_graph,
        positions=p_positions)

    # TODO: check whether sight cell incidences coincide with those of face
    # TODO: return rerun selection if needed

    return


def select_outer_face_sight_cells(outer_faces, sorted_outer_edges, is_cycle, target_vertices, graph, positions):

    # Identify all sight cells in the outer face
    outer_bounds = get_embedding_square(graph=graph, positions=positions, scaler=3)
    outer_sight_cells, outer_edge_map, outer_graph, outer_positions = get_outer_face_sight_cells(
        selected_faces=outer_faces,
        ordered_face_edges=sorted_outer_edges,
        graph=graph,
        positions=positions,
        is_cycle=is_cycle,
        bounds=outer_bounds)

    # Calculate the incidence of all sight cells to the outer face's target incident vertices
    outer_face = set().union(*outer_faces)
    outer_target_vertices = outer_face.intersection(target_vertices)
    outer_sight_cell_incidences = get_outer_face_sight_cell_incidences(sight_cells=outer_sight_cells,
                                                                       target_vertices=outer_target_vertices,
                                                                       face_edges=sorted_outer_edges,
                                                                       face_edge_map=outer_edge_map,
                                                                       positions=outer_positions)

    # Merge Outer Sight Cells with identical incidences
    outer_sight_cell_edges = get_sight_cell_edges(outer_sight_cells, outer_graph)
    merge_cells_wrapper(face_sight_cells=outer_sight_cells,
                        cells_edge_list=outer_sight_cell_edges,
                        cell_incidences=outer_sight_cell_incidences,
                        graph=outer_graph)

    # Draw Merged Embedding
    draw_graph(graph=outer_graph, positions=outer_positions)
    save_drawn_graph(f"./graph_outer_merged.png")

    # Select Sight Cells which (Together)
    selected_outer_sight_cells = select_sight_cells(cell_incidences=outer_sight_cell_incidences,
                                                    target_vertices=outer_target_vertices)
    print(f"selected cells: {selected_outer_sight_cells}")

    # Return the selected sight cells and the created graph
    return selected_outer_sight_cells, outer_graph, outer_positions


def identify_target_vertex(graph, positions):
    """

    :param graph:
    :param positions:
    :return:
    """

    # Locate all edge crossings and rank all vertices accordingly
    edge_crossings, vertex_crossings = locate_edge_crossings(graph=graph, positions=positions)
    target_vertex = get_target_vertex(vertex_crossings=vertex_crossings, graph=graph)
    target_vertex_adjacency = list(graph.neighbors(target_vertex))

    # Remove the target vertex and identify all remaining edge crossings
    remaining_graph, remaining_positions = remove_target_vertex(graph=graph,
                                                                positions=positions,
                                                                target_vertex=target_vertex)
    remaining_edge_crossings = get_remaining_edge_crossings(graph=graph,
                                                            edge_crossings=edge_crossings,
                                                            target_vertex=target_vertex)

    return target_vertex, target_vertex_adjacency, remaining_graph, remaining_positions, remaining_edge_crossings


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Command Line Arguments
    cmd_args = sys.argv
    n_vertices, m_edges, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

    # Define Input Parameters
    embedding = "kamada_kawai"

    # Diagnostics Files
    diagnostics_directory = "./output/diagnostics"
    diagnostics_file = f"barabasi_albert_{n_vertices}_{m_edges}_{seed}"

    # TESTS ------------------------------------------------------------------------------------------------------------

    # Create Output Directory
    output_directory = "./drawings/tests/test_test"
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    # Create or Load simulated graph
    print("\nCreation and Embedding of Graph")
    # graph = create_barabasi_albert_graph(n=n_vertices, m=m_edges, seed=seed)
    # positions = embed_graph(graph=graph, embedding="kamada_kawai", n_iter=None, seed=None)
    # Specify vertices and edges
    # todo: the example below causes floating point crashes as all their x and y points are identical
    # coordinates = [(0, 0), (1, 2), (2, 0), (3, 2), (4, 0), (5, 3), (4, 1), (3, 3), (2, 1), (1, 3)]
    coordinates = [(0, 2), (1, 0), (2, 1), (3, 0), (4, 2), (2, 4)]

    vertices = range(0, len(coordinates))
    edges = ((index, (index + 1) % len(vertices)) for index in range(0, len(vertices)))

    more_coordinates = [(-2, 1.5), (-1, 1.5), (-1, 0.5)]
    more_vertices = range(len(coordinates), len(coordinates) + len(more_coordinates))

    more_edges = ((more_vertices[index], more_vertices[(index + 1) % len(more_vertices)])
                  for index in range(0, len(more_vertices)))
    custom_edges = [(2, 5), (0, 4), (2, 0)]
    target_edges = [(9, 6), (9, 4), (9, 2), (9, 5)]

    # Create Graph
    graph = nx.Graph()
    for vertex in vertices:
        graph.add_node(vertex)
    for vertex in more_vertices:
        graph.add_node(vertex)
    v_index = max(graph.nodes) + 1
    graph.add_node(v_index)
    for edge in target_edges:
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)


    for edge in edges:
        print(f"edge: {edge}")
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    for edge in more_edges:
        print(f"edge: {edge}")
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    for edge in custom_edges:
        print(f"edge: {edge}")
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    positions = {vertices[index]: np.array(coordinates[index]) for index in range(0, len(coordinates))}
    positions.update(
        {more_vertices[index]: np.array(more_coordinates[index]) for index in range(0, len(more_vertices))})
    positions.update({v_index: np.array((0.0, 1.0))})
    print(graph.edges)

    # MAIN -------------------------------------------------------------------------------------------------------------

    # Draw Initial Embedding
    draw_graph(graph=graph, positions=positions)
    save_drawn_graph(f"{output_directory}/graph_0.png")

    # Identify Target and Remove it from the embedding
    print("\nIdentify Target Vertex and Remove it from the Embedding")
    target_vertex, target_adjacency, r_graph, r_positions, r_crossings = identify_target_vertex(
        graph=graph, positions=positions)

    # Draw the remaining graph
    draw_graph(graph=r_graph, positions=r_positions)
    save_drawn_graph(f"{output_directory}/graph_1.png")

    # Planarize Graph after removal of target vertex
    print("\nPlanarize Remaining Graph after target removal")
    p_graph, p_positions, virtual_edge_set = planarize_graph(graph=r_graph,
                                                             positions=r_positions,
                                                             edge_crossings=r_crossings)

    # Draw the planarized graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"{output_directory}/graph_2.png")

    # Select the faces within which to embed the split vertices
    print("\nSelect Face")
    select_embedding_faces(p_graph=p_graph,
                           p_positions=p_positions,
                           target_vertices=target_adjacency)

    sys.exit()

    # # Locate faces and best two for target face
    # TODO: find outer face
    faces = find_all_faces(plane_graph, plane_positions)
    print(f"faces: {faces}")
    face_edge_map = build_face_to_edge_map(plane_graph, faces)
    print(f"face edge map: {face_edge_map}")
    face_incidences = find_face_vertex_incidence(faces, target_vertices)
    print(f"face incidences: {face_incidences}")
    ordered_face_edges = get_ordered_face_edges(faces, plane_graph)
    print(f"ordered_face_edges: {ordered_face_edges}")
    # Get Sight Cells FOR A SELECTION OF TWO FACES
    sight_cells, edge_map = get_face_sight_cells(selected_faces=faces,
                                                 ordered_face_edges=ordered_face_edges,
                                                 graph=plane_graph,
                                                 positions=plane_positions,
                                                 bounds=((-6, -6), (-6, 6), (6, 6), (6, -6)))
    print("---------------------------------------------------------------------")
    print(f"sight cells: {sight_cells}")
    print(f"edge map: {edge_map}")

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    save_drawn_graph(f"{output_directory}/sight_cell_line_segments_1.5.png")

    sight_cell_incidences = get_faces_sight_cell_incidences(sight_cells=sight_cells,
                                                            target_vertices=target_vertices,
                                                            face_edges=ordered_face_edges,
                                                            face_edge_map=edge_map,
                                                            positions=plane_positions)
    print("----------------------------------------------------------------------")
    print(f"sight cell incidences: {sight_cell_incidences}")

    print("-----------------------------------------------------------------------")
    # Get the edge of all identified sight cells and merge sets with identical incidences
    sight_cell_edges = get_sight_cells_edge_sets(sight_cells, plane_graph)
    print(f"sight cell edges: {sight_cell_edges}")
    merge_all_face_cells(sight_cells, sight_cell_edges, sight_cell_incidences, plane_graph)

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    save_drawn_graph(f"{output_directory}/sight_cell_line_segments_1.75.png")

    sys.exit()

    #


    #
    sight_cell_edges = get_sight_cells_edge_sets(sight_cells, plane_graph)


    print(f"\nSight Cells: {sight_cells}")
    print(f"\nSight Edges: {sight_cell_edges}")
    print(f"\nIncidences:  {sight_cell_incidences}")

    # Find best set of sight cells per face
    selected_cells = select_sight_cells(sight_cells, sight_cell_incidences)
    print(f"\nSelected:    {selected_cells}")

    # Check whether selected faces match sight cell incidence
    rerank = match_cell_and_face_incidence(face_incidences=face_incidences,
                                           selected_sight_cell_incidences=selected_cells)

    # Draw and Save Planar, Convex-Face Graph
    planar_drawing_start_time = timeit.default_timer()
    draw_graph(graph=plane_graph, positions=plane_positions)
    save_drawn_graph(f"{output_directory}/sight_cell_line_segments_2.png")
    planar_drawing_time = timeit.default_timer() - planar_drawing_start_time

    sys.exit()

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
    print(f"faces: {faces}")
    #
    face_incidences = find_face_vertex_incidence(faces, target_vertex_adjacency)
    ordered_face_edges = get_ordered_face_edges(faces, plane_graph)
    print(f"face incidences: {face_incidences}")
    print(f"ordered: {ordered_face_edges}")
    print(f"faces: {faces}")

    # Find the outer face(s)
    # TODO: redo the outer face and account for branches
    # outer_face = find_outer_face(ordered_face_edges, plane_graph)
    # ordered_outer_face_edges = sort_face_edges(outer_face)
    # ordered_outer_face_vertices = get_sorted_face_vertices(ordered_outer_face_edges, is_sorted=True)
    # print(f"ordered outer face edges: {ordered_outer_face_edges}")
    # print(f"ordered outer face vertices: {ordered_outer_face_vertices}")
    # outer_face_inner_angles = calculate_face_inner_angles(ordered_outer_face_vertices, plane_positions)

    # Get Sight Cells
    sight_cells = get_face_sight_cell(faces, ordered_face_edges, plane_graph, plane_positions)
    print(f"graph: {plane_graph}")
    print(sight_cells)

    # Draw and Save Planar, Convex-Face Graph
    planar_drawing_start_time = timeit.default_timer()
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=3.5)
    save_drawn_graph(output_path)
    planar_drawing_time = timeit.default_timer() - planar_drawing_start_time

    sys.exit()
    [print(f"{key} - {ordered_face_edges[key]}") for key in ordered_face_edges.keys()]

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
    # todo: store all line-segments between runs, update them when new split places, and delete when vertex is deleted

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
