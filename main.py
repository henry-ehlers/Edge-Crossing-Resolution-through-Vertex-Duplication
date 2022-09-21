import copy

import matplotlib.pyplot as plt
import networkx as nx
import itertools as it
import pandas as pd
import numpy as np
import timeit
import sys

from src.sight_cells import *
from src.faces import *


def select_embedding_faces(incidence_table, target_vertices):
    incidence_matrix = get_incidence_matrix(incidence_table=incidence_table,
                                            targets=target_vertices)
    print(incidence_matrix)
    # TODO: resolve ties
    selected_entries = ilp_choose_face(visibility_matrix=incidence_matrix.to_numpy(dtype=int))
    print(f"selected entry indices: {selected_entries}")
    print(incidence_matrix.iloc[selected_entries])
    return selected_entries


def get_incidence_matrix(incidence_table, targets=None):

    # Initialize the incidence matrix as an empty numpy array
    entries = incidence_table["identifier"].tolist()
    targets = targets if targets is not None else list(set(incidence_table["incidence"].tolist()))
    incidence_matrix = np.empty(shape=(incidence_table.shape[0], len(targets)), dtype=int)

    # Fill the array with 1 if the target vertex of that column is present
    for row in range(0, incidence_table.shape[0]):
        for col in range(0, len(targets)):
            incidence_matrix[row][col] = 1 if targets[col] in incidence_table.loc[row, "incidence"] else 0

    # Return the matrix
    return pd.DataFrame(incidence_matrix, columns=targets, index=entries)


def get_outer_face_sight_cells(outer_faces, sorted_outer_edges, is_cycle, target_vertices, graph, positions, bounds):

    # Identify all sight cells in the outer face
    sight_cells, edge_map, connected_vertices, o_graph, o_positions = find_outer_face_sight_cells(
        selected_faces=outer_faces,
        ordered_face_edges=sorted_outer_edges,
        graph=graph,
        positions=positions,
        is_cycle=is_cycle,
        bounds=bounds)

    # Calculate the incidence of all sight cells to the outer face's target incident vertices
    outer_face = set().union(*outer_faces)
    outer_target_vertices = outer_face.intersection(target_vertices)
    cell_incidences = get_outer_face_sight_cell_incidences(sight_cells=sight_cells,
                                                           target_vertices=outer_target_vertices,
                                                           face_edges=sorted_outer_edges,
                                                           face_edge_map=edge_map,
                                                           positions=o_positions)

    # Merge Outer Sight Cells with identical incidences and Update all data structures
    outer_sight_cell_edges = get_sight_cell_edges(sight_cells, o_graph)
    sight_cells, vertex_map = merge_cells_wrapper(face_sight_cells=sight_cells,
                                                  cell_incidences=cell_incidences,
                                                  cells_edge_map=edge_map,
                                                  cells_edge_list=outer_sight_cell_edges,
                                                  positions=o_positions,
                                                  graph=o_graph)

    # DEBUG
    print(f"\n cells MAP:")
    [print(f"{cell}") for cell in cell_incidences]
    print(f"----------------------------------------------------------")
    print(f"\n incidences MAP:")
    [print(f"{cell} - {cell_incidences[cell]}") for cell in cell_incidences.keys()]
    print(f"----------------------------------------------------------")
    print(f"\n face sight cells MAP:")
    [print(f"{key}") for key in sight_cells]
    print(f"----------------------------------------------------------")
    print(f"\n face vertex map:")
    [print(f"{key} -> {item}") for key, item in vertex_map.items()]
    print(f"----------------------------------------------------------")

    # Get Sorted Incidence Table of sight cells and their incidences
    cell_incidence_table = get_incidence_table(incidences=cell_incidences,
                                               entry_type="cell",
                                               outer=True)

    # Return Everything if no single sight cell can realize the incidence of the face
    new_graph_object = {"cells":           sight_cells,
                        "incidences":      cell_incidences,
                        "connected_nodes": connected_vertices,
                        "edge_map":        edge_map,
                        "vertex_map":      vertex_map,
                        "graph":           o_graph,
                        "positions":       o_positions}
    return cell_incidence_table, new_graph_object


def get_incidence_table(incidences, entry_type="face", outer=False):

    # Create lists for each of the 4 columns
    entry_list = [entry for entry in incidences.keys()]
    incidence_list = [incidences[entry] for entry in entry_list]
    type_list = [entry_type] * len(entry_list)
    outer_list = [outer] * len(entry_list)

    # Create n x 4 data table, where n = the number of faces
    incidence_table = pd.DataFrame({"type":       type_list,
                                    "outer":      outer_list,
                                    "identifier": entry_list,
                                    "incidence":  incidence_list})

    # Return the pandas data table
    return incidence_table


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


def get_outer_face(sorted_inner_face_edges, graph, positions):

    # Find Outer Face
    outer_faces, outer_face_edges, is_cycle = find_outer_face(sorted_inner_face_edges, graph, positions)

    #
    # outer_face_identifier = frozenset(set.union(*[set(outer_face) for outer_face in outer_faces]))
    sorted_outer_face_edges = {outer_face: sort_face_edges(outer_face_edges[outer_face])
                               for outer_face in outer_faces if is_cycle[outer_face]}

    return outer_faces, sorted_outer_face_edges, is_cycle


def decompose_outer_face(sorted_inner_face_edges, target_vertices, graph, positions, bounds):

    # Find the Graph's Outer face(s)
    outer_faces, sorted_outer_face_edges, is_cycle = get_outer_face(
        sorted_inner_face_edges=sorted_inner_face_edges, graph=graph, positions=positions)

    # Decompose the Graph's Outer face
    cell_incidence_table, new_graph_object = get_outer_face_sight_cells(outer_faces=outer_faces,
                                                                        sorted_outer_edges=sorted_outer_face_edges,
                                                                        is_cycle=is_cycle,
                                                                        target_vertices=target_vertices,
                                                                        graph=graph,
                                                                        positions=positions,
                                                                        bounds=bounds)

    # Return the decomposed outer face subgraph
    return cell_incidence_table, new_graph_object


def update_graph_with_sight_cells(graph, positions, cell_graph, cell_positions, new_edge_map):
    print(f'\nOld Edges: {graph.edges()}')
    print(f"Old Nodes: {graph.nodes()}")

    print(f'\nNew Edges: {cell_graph.edges()}')
    print(f"New Nodes: {cell_graph.nodes()}")

    # Remove edges replaced with sets of virtual edges, should they already exist in the graph
    [graph.remove_edge(u=e[0], v=e[1]) for e in new_edge_map.keys() if graph.has_edge(u=e[0], v=e[1])]

    # Update the position list of vertices (i.e. add new ones; old ones do not change)
    positions.update(cell_positions)

    # Update the graph with the vertices and virtual edges
    graph.update(cell_graph)


def get_inner_faces(target_vertices, graph, positions):

    # Identify the graph's inner faces
    inner_faces = find_inner_faces(graph=graph,
                                   positions=positions)

    # Identify each face's incidence (based not on visibility)
    inner_faces_incidences = find_face_vertex_incidence(faces=inner_faces,
                                                        target_vertices=target_vertices)

    # Create Pandas Data Table of Face Incidences
    inner_incidence_table = get_incidence_table(incidences=inner_faces_incidences,
                                                entry_type="face",
                                                outer=False)

    # Get the Face's sorted Edges
    sorted_face_edges = get_ordered_face_edges(faces=inner_faces,
                                               graph=graph)

    # Return both the incidence table and the sorted edges
    return inner_incidence_table, sorted_face_edges


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
    coordinates = [(0.001, 2.003), (1.001, 0.005), (2.003, 1.002), (3.004, 0.002),
                   (4.003, 2.006), (2.001, 4.004)]

    vertices = range(0, len(coordinates))
    edges = ((index, (index + 1) % len(vertices)) for index in range(0, len(vertices)))

    more_coordinates = [(-2.0004, 1.5004), (-1.16, 1.55), (-1.004, 0.507)]
    more_vertices = range(len(coordinates), len(coordinates) + len(more_coordinates))

    more_edges = ((more_vertices[index], more_vertices[(index + 1) % len(more_vertices)])
                  for index in range(0, len(more_vertices)))
    custom_edges = [(2, 5), (0, 4), (2, 0)]
    target_edges = [(9, 6), (9, 4), (9, 2), (9, 5)]

    # Create Graph
    graph = nx.Graph()
    for vertex in vertices:
        graph.add_node(vertex, real=1)
    for vertex in more_vertices:
        graph.add_node(vertex, real=1)
    v_index = max(graph.nodes) + 1
    graph.add_node(v_index, real=1)
    for edge in target_edges:
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)

    for edge in edges:
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    for edge in more_edges:
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    for edge in custom_edges:
        graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    positions = {vertices[index]: np.array(coordinates[index]) for index in range(0, len(coordinates))}
    positions.update(
        {more_vertices[index]: np.array(more_coordinates[index]) for index in range(0, len(more_vertices))})
    positions.update({v_index: np.array((0.0, 1.0))})

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
    p_graph, p_positions = copy.deepcopy(r_graph), copy.deepcopy(r_positions)
    virtual_edge_set = planarize_graph(graph=p_graph,
                                       positions=p_positions,
                                       edge_crossings=r_crossings,
                                       largest_index=max(graph.nodes))

    # Draw the planarized graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"{output_directory}/graph_2.png")

    # Get Inner Faces
    print(f"\nIdentify the Inner Faces")
    inner_face_incidence, sorted_inner_face_edges = get_inner_faces(target_vertices=target_adjacency,
                                                                    graph=p_graph,
                                                                    positions=p_positions)
    inner_face_edge_map = get_ordered_face_edges(faces=inner_face_incidence['identifier'].tolist(),
                                                 graph=p_graph)
    print(f"inner face edge map: {inner_face_edge_map}")
    print(inner_face_incidence)
    print(f"\n sorted inner edges: {sorted_inner_face_edges}")

    # Decompose the outer face into sight cells and update the planar graph
    print("\nDecompose The Outer Face")
    d_graph, d_positions = copy.deepcopy(p_graph), copy.deepcopy(p_positions)
    outer_bounds = get_embedding_square(graph=p_graph,
                                        positions=p_positions,
                                        scaler=1.2)
    outer_cell_incidence, cell_graph_object = decompose_outer_face(sorted_inner_face_edges=sorted_inner_face_edges,
                                                                   graph=p_graph,
                                                                   positions=p_positions,
                                                                   target_vertices=target_adjacency,
                                                                   bounds=outer_bounds)
    update_graph_with_sight_cells(graph=d_graph,
                                  positions=d_positions,
                                  cell_graph=cell_graph_object["graph"],
                                  cell_positions=cell_graph_object["positions"],
                                  new_edge_map=cell_graph_object["edge_map"])

    draw_graph(graph=d_graph, positions=d_positions)
    save_drawn_graph(f"{output_directory}/graph_3.png")

    # Create line-segments between all vertices now already connected by edges or virtual edge sets
    print(f"\nUpdate Inner Face")
    update_faces_with_edge_map(inner_face_incidence, sorted_inner_face_edges, cell_graph_object["edge_map"])
    print(inner_face_incidence)
    print()
    print(sorted_inner_face_edges)

    # Select the targets within which to embed split vertices
    print(f"\nSelect Embedding Cells/Faces")
    incidence_table = pd.concat(objs=[inner_face_incidence, outer_cell_incidence],
                                ignore_index=True,
                                axis=0)
    selected_cells = select_embedding_faces(incidence_table, target_adjacency)
    print(f"selected cells: {selected_cells}")
    print(f"index: {selected_cells.index}")
    selected_faces = [incidence_table.at[row, "identifier"] for row in selected_cells]

    print(f"selected faces: {selected_faces}")

    print(f"\nDraw All-to-All line segments")
    [virtual_edge_set.pop(cell) for cell in list(virtual_edge_set.keys()) if not virtual_edge_set[cell]]
    edge_map = {**virtual_edge_set, **cell_graph_object["edge_map"]}
    for edge in list(edge_map.keys()):
        edge_map[frozenset(edge)] = edge_map[edge]
        edge_map.pop(edge)
    print(f"\nvirtual edge set:")
    [print(f"{edge} - {edge_map[edge]}") for edge in edge_map.keys()]
    print(f"\n already extended:")
    [print(f"{vertex}: {cell_graph_object['connected_nodes'][vertex]}") for vertex in cell_graph_object['connected_nodes'].keys()]
    s_graph, s_positions, s_edge_map = draw_all_line_segments(graph=d_graph,
                                                              positions=d_positions,
                                                              virtual_edge_set=edge_map,
                                                              bounds=outer_bounds,
                                                              already_extended=cell_graph_object['connected_nodes'])

    # Draw the segment graph
    draw_graph(graph=s_graph, positions=s_positions)
    save_drawn_graph(f"{output_directory}/graph_4.png")



    print(f"\nCull Non-Selected Line Segments")
    print(f"cells: {incidence_table['identifier'].tolist()}")
    # TODO: the cells of the decomposed outer face have not been updated, since we commented out the vertex deletion
    #  subsequently, the ordered_face_edges function is looking for edges which no longer exist
    #  i.e. edges that were merge edges -> these vertices (singletons) need to be removed from the incidence set
    complete_cell_edge_map = get_ordered_face_edges(faces=incidence_table['identifier'].tolist(),
                                                    graph=d_graph)
    [print(f"{cell} - {complete_cell_edge_map[cell]}") for cell in complete_cell_edge_map.keys()]

    print(f"\nFace edge map:")
    print(complete_cell_edge_map)
    c_graph, c_positions, intersection_map = cull_all_line_segment_graph(
        target_faces=selected_faces,
        face_edge_map=complete_cell_edge_map,
        face_vertex_map=cell_graph_object["vertex_map"],
        segment_edge_map=s_edge_map,
        graph=s_graph,
        positions=s_positions)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{output_directory}/graph_5.png")

    subface_edge_set, subface_vertex_map = create_subface_graph(
        graph=c_graph,
        positions=c_positions,
        target_faces=selected_faces,
        face_intersection_map=intersection_map)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{output_directory}/graph_6.png")

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
