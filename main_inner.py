import copy

import matplotlib.pyplot as plt
import networkx as nx
import itertools as it
import pandas as pd
import numpy as np
import datetime
import pickle
import sys

from src.graph_simulation import *
from src.vertex_splitting import *
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


def weighted_select_embedding_faces(incidence_table, edge_length_table, target_vertices):
    print(incidence_table)
    print(edge_length_table)
    # input("HERE")
    incidence_matrix = get_incidence_matrix(incidence_table=incidence_table,
                                            targets=target_vertices)
    print(incidence_matrix)
    edge_length_matrix = get_edge_length_matrix(edge_length_table, incidence_table, target_vertices)
    print(edge_length_matrix)
    # input("CHJELJK")
    assert (incidence_matrix.shape[0] == edge_length_matrix.shape[0]) and \
           (incidence_matrix.shape[1] == edge_length_matrix.shape[1]), \
        "Error in face selection. Incidence Matrix and Edge Length Matrix must of of same dimensions."

    print(incidence_matrix)
    # TODO: resolve ties
    selected_entries = ilp_choose_weighted_face(visibility_matrix=incidence_matrix.to_numpy(dtype=int),
                                                edge_length_dif=edge_length_matrix)
    print(f"selected entry indices: {selected_entries}")
    print(incidence_matrix.iloc[selected_entries])
    # input("selected faces")
    return selected_entries


def calculate_edge_length_weights(sight_cells, targets, positions, target_length):
    edge_lengths = [[None] * len(sight_cells) for target in targets]
    print(len(edge_lengths))
    print(len(edge_lengths[0]))

    for cell_index, sight_cell in enumerate(sight_cells):
        print(f"cell {sight_cell}")
        face_vertex_positions = [positions[v] for v in list(sight_cell)]
        print(f"positions: {face_vertex_positions}")
        centroid = calculate_face_centroid(face_vertex_positions)
        print(f"centroid: {centroid}")
        for target_index, target in enumerate(targets):
            target_position = positions[target]
            print(f"target: {target} @ {target_position}")
            edge_lengths[target_index][cell_index] = abs(squared_distance(centroid, target_position) - target_length)
    [print(edge_lengths[i]) for i in range(0, len(edge_lengths))]

    edge_length_table = pd.DataFrame({t: edge_lengths[t_index] for t_index, t in enumerate(targets)})
    edge_length_table["identifier"] = sight_cells

    return edge_length_table


def get_face_centroids(faces: [set], positions):
    centroids = {face: None for face in faces}
    for face in faces:
        face_vertex_positions = [positions[v] for v in list(face)]
        print(f"positions: {face_vertex_positions}")
        centroids[face] = calculate_face_centroid(face_vertex_positions)
    return centroids


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


def get_edge_length_matrix(edge_length_table, incidence_table, targets):
    edge_length_matrix = np.empty(shape=(incidence_table.shape[0], len(targets)), dtype=float)
    print(edge_length_table[targets[0]])
    sight_cells = incidence_table["identifier"].tolist()
    for cell_index, cell in enumerate(sight_cells):
        print(f"Cell {cell} @ {cell_index}")
        for target_index, target in enumerate(targets):
            match_index = edge_length_table.index[edge_length_table["identifier"] == cell].tolist()[0]
            print(f"Target {target} @ {target_index}")
            print(edge_length_table["identifier"] == cell)
            print(f'index: {match_index}')
            edge_length = edge_length_table.at[match_index, target]
            edge_length_matrix[cell_index, target_index] = edge_length

    return edge_length_matrix


def get_outer_face_sight_cells(outer_faces, sorted_outer_edges, is_cycle, target_vertices, graph, positions, bounds):

    # Identify all sight cells in the outer face
    sight_cells, ordered_cell_edges, edge_map, connected_vertices, o_graph, o_positions = find_outer_face_sight_cells(
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
    sight_cells, ordered_cell_edges, vertex_map = merge_cells_wrapper(face_sight_cells=sight_cells,
                                                                      cell_incidences=cell_incidences,
                                                                      cells_edge_map=edge_map,
                                                                      ordered_cell_edges=ordered_cell_edges,
                                                                      positions=o_positions,
                                                                      graph=o_graph)

    # Get Sorted Incidence Table of sight cells and their incidences
    cell_incidence_table = get_incidence_table(incidences=cell_incidences,
                                               entry_type="cell",
                                               outer=True)

    # Return Everything if no single sight cell can realize the incidence of the face
    new_graph_object = {"cells": sight_cells,
                        "incidences": cell_incidences,
                        "connected_nodes": connected_vertices,
                        "ordered_cycle_edges": ordered_cell_edges,
                        "edge_map": edge_map,
                        "vertex_map": vertex_map,
                        "graph": o_graph,
                        "positions": o_positions}

    return cell_incidence_table, new_graph_object


def get_incidence_table(incidences, entry_type="face", outer=False):
    # Create lists for each of the 4 columns
    entry_list = [entry for entry in incidences.keys()]
    incidence_list = [incidences[entry] for entry in entry_list]
    type_list = [entry_type] * len(entry_list)
    outer_list = [outer] * len(entry_list)

    # Create n x 4 data table, where n = the number of faces
    incidence_table = pd.DataFrame({"type": type_list,
                                    "outer": outer_list,
                                    "identifier": entry_list,
                                    "incidence": incidence_list})

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
    [print(f"{key} - {item}") for key, item in edge_crossings.items()]
    print()
    print(vertex_crossings)
    if len(edge_crossings) == 0:
        return None, None, None, None, None

    target_vertex = get_target_vertex(vertex_crossings=vertex_crossings, graph=graph)
    print(f"target vertex: {target_vertex}")
    # input("---------")
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
    print(f"\nOuter Face: {outer_faces}")
    print(f"\nOuter Face Edges: {outer_face_edges}")
    #
    # outer_face_identifier = frozenset(set.union(*[set(outer_face) for outer_face in outer_faces]))
    sorted_outer_face_edges = {outer_face: sort_face_edges(outer_face_edges[outer_face])
                               for outer_face in outer_faces if is_cycle[outer_face]}

    return outer_faces, sorted_outer_face_edges, is_cycle


def decompose_outer_face(sorted_inner_face_edges, target_vertices, graph, positions, bounds):
    # Find the Graph's Outer face(s)
    outer_faces, sorted_outer_face_edges, is_cycle = get_outer_face(
        sorted_inner_face_edges=sorted_inner_face_edges,
        graph=graph,
        positions=positions)
    print(f'outer faces: {outer_faces}')
    # input("......")

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
    for virtualized_edge in new_edge_map.keys():
        virtualized_edge = tuple(virtualized_edge)
        if graph.has_edge(u=virtualized_edge[0], v=virtualized_edge[1]):
            graph.remove_edge(u=virtualized_edge[0], v=virtualized_edge[1])

    # Update the position list of vertices (i.e. add new ones; old ones do not change)
    positions.update(cell_positions)

    # Update the graph with the vertices and virtual edges
    graph.update(cell_graph)


def get_inner_faces(target_vertices: [int], graph, positions, outer_bounds):

    # Identify the graph's inner faces
    # input("find_inner_faces")
    inner_faces, sorted_inner_face_vertices, sorted_inner_face_edges = find_inner_faces(graph=graph,
                                                                                        positions=positions)
    # print(f"\nINNER FACES:")
    [print(face) for face in inner_faces]

    # Decompose Inner Faces into Sight Cells by sight line extension
    # input("find_inner_face_sight_cells")
    cells, ordered_cell_vertices, ordered_cell_edges, face_edge_map, connected_vertex_map = find_inner_face_sight_cells(
        inner_faces=inner_faces,
        ordered_face_edges=sorted_inner_face_edges,
        graph=graph,
        positions=positions,
        bounds=outer_bounds)

    print(f"\n connected nodes:")
    [print(f"cell: {cell} - {vertices}") for cell, vertices in connected_vertex_map.items()]

    print(f"\nINNER CELLS:")
    [print(cell) for cell in cells]
    # input(f"inner CELLS")

    print(f"ordered cell vertices:")
    [print(f"cell: {cell} - {vertices}") for cell, vertices in ordered_cell_vertices.items()]
    # input("...")

    print(f"ordered cell edges:")
    [print(f"cell: {cell} - {edges}") for cell, edges in ordered_cell_edges.items()]
    # input("...")

    print(f"face edge map:")
    [print(f"cell: {cell} - {edges}") for cell, edges in face_edge_map.items()]

    # Create Pandas Data Table of Face Incidences
    # input("find_face_vertex_incidence")
    convex_faces = inner_faces.intersection(cells)
    inner_faces_incidences = find_face_vertex_incidence(faces=convex_faces,
                                                        target_vertices=target_vertices)

    # Get Incidence of Sight Cells Identified
    actual_cells = cells - convex_faces
    print(f"actual cells: {actual_cells}")
    # input("get_inner_face_sight_cell_incidences")
    # cell_edges = get_sight_cell_edges(actual_cells, graph)
    inner_cells_incidence = get_inner_face_sight_cell_incidences(sight_cells=actual_cells,
                                                                 target_vertices=target_vertices,
                                                                 face_edges=sorted_inner_face_edges,
                                                                 face_edge_map=face_edge_map,
                                                                 positions=positions)

    print(f"\ninner_cells_incidence")
    print(inner_cells_incidence)
    # input("merge_cells_wrapper")
    sight_cells, ordered_cell_edges, vertex_map = merge_cells_wrapper(face_sight_cells=actual_cells,
                                                                      cell_incidences=inner_cells_incidence,
                                                                      cells_edge_map=face_edge_map,
                                                                      ordered_cell_edges=ordered_cell_edges,
                                                                      positions=positions,
                                                                      graph=graph)
    print(f"sight_cells: {sight_cells}")
    print(f"ordered edges:")
    ordered_edges = {**{face: sorted_inner_face_edges[face] for face in convex_faces},
                     **{cell: ordered_cell_edges[cell] for cell in ordered_cell_edges.keys()}}
    [print(f"{cell} - {item}") for cell, item in ordered_edges.items()]
    print(f"\nedge map:")
    [print(f"{cell} - {item}") for cell, item in face_edge_map.items()]
    print(f"\nvertex map:")
    [print(f"{cell} - {item}") for cell, item in vertex_map.items()]
    print(f"ordered face edges:")
    [print(f"cell: {face} - {edges}") for face, edges in sorted_inner_face_edges.items()]
    print(f"ordered cell edges:")
    [print(f"cell: {cell} - {edges}") for cell, edges in ordered_cell_edges.items()]
    print(f"ordered cell edges:")
    [print(f"cell: {key} - {edges}") for key, edges in ordered_edges.items()]
    # input("...")

    # Get Combined Incidence Table
    cell_incidence_table = get_incidence_table(incidences=inner_cells_incidence, entry_type="cell", outer=False)
    print(cell_incidence_table)
    # input("...")
    face_incidence_table = get_incidence_table(incidences=inner_faces_incidences, entry_type="face", outer=False)
    print(face_incidence_table)
    # input("...")
    inner_incidence_table = pd.concat(objs=[face_incidence_table, cell_incidence_table], ignore_index=True)
    print(f"\ninner_incidence_table:")
    print(inner_incidence_table)

    print(f"\nface edge map:")
    print(face_edge_map)
    vertex_map.update({face: face for face in convex_faces})
    print(f"\nconnected_vertex_map:")
    print(connected_vertex_map)

    # Return Everything if no single sight cell can realize the incidence of the face
    new_graph_object = {"cells": sight_cells,
                        "incidences": inner_cells_incidence,
                        "connected_nodes": connected_vertex_map,
                        "ordered_cycle_edges": ordered_edges,
                        "edge_map": face_edge_map,
                        "vertex_map": vertex_map}

    # Return both the incidence table and the sorted edges
    return inner_incidence_table, new_graph_object


def split_vertex(graph, positions, labels, target_edge_length=1.0, drawing_directory="."):
    # Draw Initial Embedding
    draw_graph(graph=graph, positions=positions, labels=labels)
    save_drawn_graph(f"{drawing_directory}/graph_0.png")

    # Identify Target and Remove it from the embedding
    print("\nIdentify Target Vertex and Remove it from the Embedding")
    target_vertex, target_adjacency, r_graph, r_positions, r_crossings = identify_target_vertex(
        graph=graph, positions=positions)

    print(f"target vertex: {target_vertex}")
    print(f"target adjacency: {target_adjacency}")
    # input("FLKSDFJ")
    if target_vertex is None:
        return False, graph, positions, labels

    # Draw the remaining graph
    draw_graph(graph=r_graph, positions=r_positions)
    save_drawn_graph(f"{drawing_directory}/graph_1.png")

    # Planarize Graph after removal of target vertex
    print("\nPlanarize Remaining Graph after target removal")
    rr_graph, rr_positions = copy.deepcopy(r_graph), copy.deepcopy(r_positions)
    virtual_edge_map = planarize_graph(graph=rr_graph,
                                       positions=rr_positions,
                                       edge_crossings=r_crossings,
                                       largest_index=max(graph.nodes))

    #
    p_graph, p_positions = copy.deepcopy(rr_graph), copy.deepcopy(rr_positions)
    connect_singleton_vertex_edges(p_graph, p_positions)
    # input("CONNECTING SINGLETONS")

    # Draw the remaining graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"{drawing_directory}/graph_1.5.png")

    # virtual_edge_map = planarize_graph(graph=p_graph,
    #                                    positions=p_positions,
    #                                    edge_crossings=r_crossings,
    #                                    largest_index=max(graph.nodes))
    print(f"virtual edge set:")
    [print(f"{key} - {item}") for key, item in virtual_edge_map.items()]

    # Draw the planarized graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"{drawing_directory}/graph_2.png")

    # Get Inner Faces
    print(f"\nIdentify the Inner Faces")
    outer_bounds = get_embedding_square(graph=p_graph, positions=p_positions, scaler=1.5)
    print(f"outer bounds: {outer_bounds}")
    # input("as;lkdfjals;kdf")
    inner_face_incidence, inner_graph_object = get_inner_faces(target_vertices=target_adjacency,
                                                               graph=p_graph,
                                                               positions=p_positions,
                                                               outer_bounds=outer_bounds)
    print("POSITIONS:")
    [print(f"{key} - {value}") for key, value in p_positions.items()]
    print("INCIDENCES")
    print(inner_face_incidence)
    # input("")

    # Draw the inner sight cell decomposed graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"{drawing_directory}/graph_2.png")

    # sys.exit()

    # Decompose the outer face into sight cells and update the planar graph
    # print("\nDecompose The Outer Face")
    # outer_cell_incidence, cell_graph_object = decompose_outer_face(
    #     sorted_inner_face_edges=inner_graph_object["ordered_cycle_edges"],
    #     graph=p_graph,
    #     positions=p_positions,
    #     target_vertices=target_adjacency,
    #     bounds=outer_bounds)

    # draw_graph(graph=cell_graph_object["graph"], positions=cell_graph_object["positions"])
    # save_drawn_graph(f"{drawing_directory}/graph_2.5.png")
    #
    # d_graph, d_positions = copy.deepcopy(p_graph), copy.deepcopy(p_positions)
    # update_graph_with_sight_cells(graph=d_graph,
    #                               positions=d_positions,
    #                               cell_graph=cell_graph_object["graph"],
    #                               cell_positions=cell_graph_object["positions"],
    #                               new_edge_map=cell_graph_object["edge_map"])
    #
    # draw_graph(graph=d_graph, positions=d_positions)
    # save_drawn_graph(f"{drawing_directory}/graph_3.png")

    # TODO: fix the projection to now work with the data structures of edges
    # TODO: fix subface identification to work with new face detection function and produced data structure

    # -------------------------------------------------------------------------------------------------------
    # input("\nCONNECTED NODES ----------------------------------------------------------------------------")
    # print(f"\nA")
    # [print(f"{key} - {items}") for key, items in inner_graph_object['connected_nodes'].items()]
    # print(f"\nB")
    # [print(f"{key} - {items}") for key, items in cell_graph_object['connected_nodes'].items()]
    # connected_nodes = merge_connected_nodes([inner_graph_object['connected_nodes'],
    #                                          cell_graph_object['connected_nodes']])
    # print(f"\nC")
    # [print(f"{key} - {items}") for key, items in connected_nodes.items()]

    # -------------------------------------------------------------------------------------------------------
    # input("\nEDGE MAP ----------------------------------------------------------------------------------")
    # print(f"\nA")
    # [print(f"{key} - {items}") for key, items in inner_graph_object['edge_map'].items()]
    # print(f"\nB")
    # [print(f"{key} - {items}") for key, items in cell_graph_object['edge_map'].items()]
    # complete_edge_map = merge_edge_map(old_edge_map=inner_graph_object['edge_map'],
    #                                    new_edge_map=cell_graph_object['edge_map'])
    # print(f"\nC")
    # [print(f"{key} - {items}") for key, items in complete_edge_map.items()]

    # input("CHECK")

    # -------------------------------------------------------------------------------------------------------
    # print("\nINCIDENCE TABLE ---------------------------------------------------------------------------")
    #
    # # Create line-segments between all vertices now already connected by edges or virtual edge sets
    # print(f"\nUpdate Inner Face")
    # sorted_inner_face_edges = inner_graph_object["ordered_cycle_edges"]
    # update_faces_with_edge_map(face_incidence_table=inner_face_incidence,
    #                            sorted_face_edges=sorted_inner_face_edges,
    #                            edge_map=complete_edge_map)
    # # input(f"\n POST UPDATE")
    #
    # # Select the targets within which to embed split vertices
    # print(f"\nSelect Embedding Cells/Faces")
    # incidence_table = pd.concat(objs=[inner_face_incidence, outer_cell_incidence],
    #                             ignore_index=True,
    #                             axis=0)
    #
    # # TODO: calculcate lengths
    # edge_length_table = calculate_edge_length_weights(sight_cells=incidence_table["identifier"].tolist(),
    #                                                   targets=target_adjacency,
    #                                                   positions=d_positions,
    #                                                   target_length=target_edge_length)
    #
    # print(incidence_table)
    # print(f"")
    # selected_cells = weighted_select_embedding_faces(incidence_table=incidence_table,
    #                                                  edge_length_table=edge_length_table,
    #                                                  target_vertices=target_adjacency)
    # # selected_cells = select_embedding_faces(incidence_table=incidence_table,
    # #                                         target_vertices=target_adjacency)
    # selected_faces = [incidence_table.at[row, "identifier"] for row in selected_cells]
    # print(f"\nselected faces: {selected_cells}")
    # print(f"\nselected faces: {selected_faces}")
    # # input(";alsdkjf;lsadkfa;sdlkfj")
    # # input("Start All Line Segmentation")

    # All-to-All Line Segments ---------------------------------------------------------------------------------

    d_graph, d_positions = p_graph, p_positions

    connected_nodes = inner_graph_object['connected_nodes']
    complete_edge_map = inner_graph_object['edge_map']
    incidence_table = inner_face_incidence

    edge_length_table = calculate_edge_length_weights(sight_cells=incidence_table["identifier"].tolist(),
                                                      targets=target_adjacency,
                                                      positions=d_positions,
                                                      target_length=target_edge_length)

    selected_cells = weighted_select_embedding_faces(incidence_table=incidence_table,
                                                     edge_length_table=edge_length_table,
                                                     target_vertices=target_adjacency)

    selected_faces = [incidence_table.at[row, "identifier"] for row in selected_cells]

    s_graph, s_positions, s_edge_map = draw_all_line_segments(graph=d_graph,
                                                              positions=d_positions,
                                                              virtual_edge_set=complete_edge_map,
                                                              bounds=outer_bounds,
                                                              already_extended=connected_nodes)
    # print(f"\nS edge Map:")
    [print(f"{k} - {v}") for k, v in s_edge_map.items()]

    # Draw the segment graph
    draw_graph(graph=s_graph, positions=s_positions)
    save_drawn_graph(f"{drawing_directory}/graph_4.png")

    # print(f"\nCull Non-Selected Line Segments")
    # print(f"cells: {incidence_table['identifier'].tolist()}")
    complete_vertex_map = inner_graph_object["vertex_map"]
    complete_face_edges = inner_graph_object["ordered_cycle_edges"]

    # print(f"A ------------------------------------------------------------------------")
    # [print(f"{key} - {item}") for key, item in inner_graph_object["vertex_map"].items()]
    # print(f"B ------------------------------------------------------------------------")
    # [print(f"{key} - {item}") for key, item in cell_graph_object["vertex_map"].items()]
    # complete_vertex_map = {**inner_graph_object["vertex_map"],
    #                        **cell_graph_object["vertex_map"]}
    # print(f"C ------------------------------------------------------------------------")
    # [print(f"{key} - {item}") for key, item in complete_vertex_map.items()]
    # input("check complete")

    # TODO: this merger MAY be incomplete? -> if an inner face's outer edge was bisected by an outer face's sight
    #  extension, this may not be accounted for
    # complete_face_edges = {**inner_graph_object["ordered_cycle_edges"],
    #                        **cell_graph_object["ordered_cycle_edges"]}
    # print(f"A ------------------------------------------------------------------------")
    # [print(f"{key} - {item}") for key, item in inner_graph_object["ordered_cycle_edges"].items()]
    # print(f"B ------------------------------------------------------------------------")
    # [print(f"{key} - {item}") for key, item in cell_graph_object["ordered_cycle_edges"].items()]
    # print(f"C ------------------------------------------------------------------------")
    # [print(f"{key} - {item}") for key, item in complete_face_edges.items()]
    # input("CHECK")

    # Cull all segments which do not intersect the two selected faces
    c_graph, c_positions, intersection_map = cull_all_line_segment_graph(
        target_faces=selected_faces,
        face_edge_map=complete_face_edges,
        face_vertex_map=complete_vertex_map,
        segment_edge_map=s_edge_map,
        graph=s_graph,
        positions=s_positions)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{drawing_directory}/graph_5.png")

    # Create Subface Graph ---------------------------------------------------------------------------------------------

    #
    subface_edge_set, subface_vertex_map = create_subface_graph(graph=c_graph,
                                                                positions=c_positions,
                                                                target_faces=selected_faces,
                                                                face_vertex_map=complete_vertex_map,
                                                                face_intersection_map=intersection_map)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{drawing_directory}/graph_6.png")

    # Identify all edge crossings in the faces, planarize the graph, and update the face's vertex sets
    print(f"\nSUBFACE CREATION")
    face_edge_crossings, face_vertex_crossings = locate_edge_crossings(graph=c_graph,
                                                                       positions=c_positions)
    plane_face_virtual_edge_map = planarize_graph(graph=c_graph,
                                                  positions=c_positions,
                                                  edge_crossings=face_edge_crossings)
    update_face_vertex_map(vertex_map=subface_vertex_map,
                           virtual_edge_map=plane_face_virtual_edge_map)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{drawing_directory}/graph_7.png")

    # Select Subfaces
    print(f"\nSUBFACE SELECTION")
    plane_graph_sub_faces = find_all_subfaces(target_faces=selected_faces,
                                              face_vertex_map=subface_vertex_map,
                                              graph=c_graph,
                                              positions=c_positions)

    # Calculate each subface's centroid
    subface_centroids = get_split_vertex_locations(positions=c_positions,
                                                   target_face_subfaces=plane_graph_sub_faces)

    # TODO: get distances to centroids
    face_centroids = get_face_centroids(selected_faces, c_positions)
    subface_distances = get_subface_distance_to_face_centroids(target_faces=selected_faces,
                                                               face_centroids=face_centroids,
                                                               subfaces=plane_graph_sub_faces,
                                                               subface_centroids=subface_centroids)
    [print(val) for key, val in subface_distances.items()]
    # input("CLKSL:KDJF:LKJ")

    # Calculate the number of edge crossing induced by connected each subface to all target neighbors
    induced_edge_crossings = calculate_induced_edge_crossings(graph=r_graph,
                                                              positions=r_positions,
                                                              centroids=subface_centroids,
                                                              target_neighbors=target_adjacency)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{drawing_directory}/graph_8.png")

    # Place Copies
    print(f"\nPLACING SPLIT VERTICES")

    # Get the number of induced edge crossings in the form of a dictionary of pandas dataframes
    induced_edge_crossing_table = get_edge_crossing_table(induced_edge_crossings=induced_edge_crossings,
                                                          target_neighbors=target_adjacency)

    #
    selected_sub_faces = select_weighted_sub_faces(sub_face_tables=induced_edge_crossing_table,
                                                   subface_distances=subface_distances,
                                                   target_faces=selected_faces,
                                                   target_vertices=target_adjacency)
    # input()
    print(f"selected sub_faces: {selected_sub_faces}")

    n_graph, n_positions = place_split_vertices(faces=selected_faces,
                                                selected_sub_faces=selected_sub_faces,
                                                centroids=subface_centroids,
                                                target_vertex=target_vertex,
                                                graph=r_graph,
                                                positions=r_positions)

    new_vertices = [v for v in n_graph.nodes() if v not in r_graph.nodes]
    print(f"new vertices: {new_vertices}")
    [labels.update({v: labels[target_vertex]}) for v in new_vertices]
    print(labels)
    print(f"labels:")
    [print(f"{vertex} - {labels}") for vertex, label in labels.items()]

    # Draw the segment graph
    draw_graph(graph=n_graph, positions=n_positions, labels=labels)
    save_drawn_graph(f"{drawing_directory}/graph_9.png")

    return True, n_graph, n_positions, labels


def get_target_edge_length(graph, positions):
    distances = [None] * len(graph.edges)
    for edge_index, edge in enumerate(list(graph.edges)):
        print(f"edge: {edge}")
        vertex_a, vertex_b = edge
        print(f"vertices: {vertex_a} and {vertex_b}")
        distances[edge_index] = squared_distance(positions[vertex_a], positions[vertex_b])
        print(f"distance: {distances[edge_index]}")
    return sum(distances) / len(distances)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Command Line Arguments
    cmd_args = sys.argv
    n_vertices, m_edges, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

    # Set Recursion limit
    # sys.setrecursionlimit(1500)

    # Define Input Parameters
    embedding = "kamada_kawai"

    # Diagnostics Files
    diagnostics_directory = "./output/diagnostics"
    simulation_type = "watts_strogatz"
    diagnostics_file = f"{simulation_type}_{n_vertices}_{m_edges}_{seed}"

    # Create Output Directory
    output_directory = f"./drawings/kamada_kawai/{simulation_type}_{n_vertices}_{m_edges}_{seed}"

    # Create or Load simulated graph
    print("\nCreation and Embedding of Graph")

    # Simulate and Embed the input graph
    graph = create_barabasi_albert_graph(n=n_vertices, m=m_edges, seed=seed, type=simulation_type)
    positions = embed_graph(graph=graph, embedding="kamada_kawai", n_iter=None, seed=None)
    labels = {node: node for node in graph.nodes}

    # Save Initial Embedding
    pickle.dump(graph, open(f"{output_directory}/graph.txt", 'wb'))
    pickle.dump(positions, open(f"{output_directory}/positions.txt", 'wb'))
    pickle.dump(labels, open(f"{output_directory}/labels.txt", 'wb'))

    # Calculate the target edge length (= the average edge length in the original embedding)
    target_edge_length = get_target_edge_length(graph, positions)

    # Calculate the number of edge crossing in the initial embedding
    initial_edge_crossings, initial_vertex_crossings = locate_edge_crossings(graph=graph, positions=positions)
    remaining_edge_crossings, remaining_vertex_crossings = initial_edge_crossings, initial_vertex_crossings
    split = True

    # Resolve at least 50% of all edge crossings of the original embedding
    iteration_number = -1
    while len(remaining_edge_crossings) > (0.5 * len(initial_edge_crossings)) and split:
        iteration_number += 1

        # Create Output Directory
        drawing_directory = f"./{output_directory}/{iteration_number}/"
        Path(drawing_directory).mkdir(parents=True, exist_ok=True)

        # Resolve edge crossings by splitting one vertex
        start_time = datetime.datetime.now()
        split, graph, positions, labels = split_vertex(graph=graph,
                                                       positions=positions,
                                                       labels=labels,
                                                       target_edge_length=target_edge_length,
                                                       drawing_directory=drawing_directory)
        time_taken = datetime.datetime.now() - start_time

        # Save graph objects and time taken
        pickle.dump(graph, open(f"{drawing_directory}/graph.txt", 'wb'))
        pickle.dump(positions, open(f"{drawing_directory}/positions.txt", 'wb'))
        pickle.dump(labels, open(f"{drawing_directory}/labels.txt", 'wb'))
        pickle.dump(time_taken, open(f"{drawing_directory}/time.txt", 'wb'))

        # Calculate the number of remaining edge crossings
        remaining_edge_crossings, remaining_vertex_crossings = locate_edge_crossings(graph=graph, positions=positions)
        # input(f"Did it split -> {split}?")
