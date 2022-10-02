import copy

import matplotlib.pyplot as plt
import networkx as nx
import itertools as it
import pandas as pd
import numpy as np
import timeit
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
    print(sight_cells)
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
    print(f"\ninner faces:")
    print(inner_faces)

    # Decompose Inner Faces
    ordered_face_edges = get_ordered_face_edges(faces=inner_faces, graph=graph)
    cells, face_edge_map, connected_vertex_map = find_inner_face_sight_cells(inner_faces=inner_faces,
                                                                             ordered_face_edges=ordered_face_edges,
                                                                             graph=graph,
                                                                             positions=positions)

    # Draw the planarized graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"./testing_inner_sight_cells.png")



    # Create Pandas Data Table of Face Incidences
    convex_faces = inner_faces.intersection(cells)
    inner_faces_incidences = find_face_vertex_incidence(faces=convex_faces, target_vertices=target_vertices)
    print(f"\ninner face incidence:")
    print(inner_faces_incidences)

    # Get Incidence of Sight Cells Identified
    actual_cells = cells - convex_faces
    cell_edges = get_sight_cell_edges(actual_cells, graph)
    inner_cells_incidence = get_inner_face_sight_cell_incidences(sight_cells=actual_cells,
                                                                 target_vertices=target_vertices,
                                                                 face_edges=ordered_face_edges,
                                                                 face_edge_map=face_edge_map,
                                                                 positions=positions)
    print(f"\ninner_cells_incidence")
    print(inner_cells_incidence)
    sight_cells, vertex_map = merge_cells_wrapper(face_sight_cells=actual_cells,
                                                  cell_incidences=inner_cells_incidence,
                                                  cells_edge_map=face_edge_map,
                                                  cells_edge_list=cell_edges,
                                                  positions=positions,
                                                  graph=graph)
    print(f"sight_cells: {sight_cells}")
    print(f"vertex map: {vertex_map}")

    # Draw the planarized graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"./testing_inner_sight_cells_merged.png")

    # Get Combined Incidence Table
    cell_incidence_table = get_incidence_table(incidences=inner_cells_incidence, entry_type="cell", outer=False)
    face_incidence_table = get_incidence_table(incidences=inner_faces_incidences, entry_type="face", outer=False)
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
                        "edge_map": face_edge_map,
                        "vertex_map": vertex_map}

    # Return both the incidence table and the sorted edges
    return inner_incidence_table, new_graph_object


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
    output_directory = f"./drawings/kamada_kawai/barabasi_albert_{n_vertices}_{m_edges}_{seed}"
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    # Create or Load simulated graph
    print("\nCreation and Embedding of Graph")
    #
    # Specify vertices and edges
    # todo: the example below causes floating point crashes as all their x and y points are identical
    # coordinates = [(0, 0), (1, 2), (2, 0), (3, 2), (4, 0), (5, 3), (4, 1), (3, 3), (2, 1), (1, 3)]
    # coordinates = [(0.001, 2.003), (1.001, 0.005), (2.003, 1.002), (3.004, 0.002),
    #                (4.003, 2.006), (2.001, 4.004)]
    #
    # vertices = range(0, len(coordinates))
    # edges = ((index, (index + 1) % len(vertices)) for index in range(0, len(vertices)))
    #
    # more_coordinates = [(-2.0004, 1.5004), (-1.16, 1.55), (-1.004, 0.507)]
    # more_vertices = range(len(coordinates), len(coordinates) + len(more_coordinates))
    #
    # more_edges = ((more_vertices[index], more_vertices[(index + 1) % len(more_vertices)])
    #               for index in range(0, len(more_vertices)))
    # custom_edges = [(2, 5), (0, 4), (2, 0)]
    # target_edges = [(9, 6), (9, 4), (9, 2), (9, 5)]
    #
    # # Create Graph
    # graph = nx.Graph()
    # for vertex in vertices:
    #     graph.add_node(vertex, real=1)
    # for vertex in more_vertices:
    #     graph.add_node(vertex, real=1)
    # v_index = max(graph.nodes) + 1
    # graph.add_node(v_index, real=1)
    # for edge in target_edges:
    #     graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    #
    # for edge in edges:
    #     graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    # for edge in more_edges:
    #     graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    # for edge in custom_edges:
    #     graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
    # positions = {vertices[index]: np.array(coordinates[index]) for index in range(0, len(coordinates))}
    # positions.update(
    #     {more_vertices[index]: np.array(more_coordinates[index]) for index in range(0, len(more_vertices))})
    # positions.update({v_index: np.array((0.0, 1.0))})
    #
    # graph.add_node(10, real=1)
    # graph.add_edge(u_of_edge=9, v_of_edge=10, real=1)
    # positions[10] = (2.49, 2.55)

    graph = create_barabasi_albert_graph(n=n_vertices, m=m_edges, seed=seed)
    positions = embed_graph(graph=graph, embedding="kamada_kawai", n_iter=None, seed=None)

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
    inner_face_incidence, inner_graph_object = get_inner_faces(target_vertices=target_adjacency,
                                                               graph=p_graph,
                                                               positions=p_positions)

    print(inner_face_incidence)

    # Get the Face's sorted Edges
    sorted_inner_face_edges = get_ordered_face_edges(faces=inner_face_incidence["identifier"].tolist(), graph=p_graph)
    print(f"\n sorted inner edges: {sorted_inner_face_edges}")

    # Draw the inner sight cell decomposed graph
    draw_graph(graph=p_graph, positions=p_positions)
    save_drawn_graph(f"{output_directory}/graph_2.png")

    # Decompose the outer face into sight cells and update the planar graph
    print("\nDecompose The Outer Face")
    outer_bounds = get_embedding_square(graph=p_graph,positions=p_positions,scaler=1.5)
    outer_cell_incidence, cell_graph_object = decompose_outer_face(sorted_inner_face_edges=sorted_inner_face_edges,
                                                                   graph=p_graph,
                                                                   positions=p_positions,
                                                                   target_vertices=target_adjacency,
                                                                   bounds=outer_bounds)

    draw_graph(graph=cell_graph_object["graph"], positions=cell_graph_object["positions"])
    save_drawn_graph(f"{output_directory}/graph_2.5.png")

    d_graph, d_positions = copy.deepcopy(p_graph), copy.deepcopy(p_positions)
    update_graph_with_sight_cells(graph=d_graph,
                                  positions=d_positions,
                                  cell_graph=cell_graph_object["graph"],
                                  cell_positions=cell_graph_object["positions"],
                                  new_edge_map=cell_graph_object["edge_map"])

    draw_graph(graph=d_graph, positions=d_positions)
    save_drawn_graph(f"{output_directory}/graph_3.png")

    # Create line-segments between all vertices now already connected by edges or virtual edge sets
    print(f"\nUpdate Inner Face")
    update_faces_with_edge_map(inner_face_incidence,
                               sorted_inner_face_edges,
                               cell_graph_object["edge_map"])
    print(f"\ninner face incidence:")
    print(inner_face_incidence)

    print(f"\nouter face incidence:")
    print(outer_cell_incidence)

    # Select the targets within which to embed split vertices
    print(f"\nSelect Embedding Cells/Faces")
    incidence_table = pd.concat(objs=[inner_face_incidence, outer_cell_incidence],
                                ignore_index=True,
                                axis=0)
    print(f"incidence table:")
    print(incidence_table)
    input("Please ENTER...")
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
    [print(f"{vertex}: {cell_graph_object['connected_nodes'][vertex]}")
     for vertex in cell_graph_object['connected_nodes'].keys()]

    connected_nodes = {**inner_graph_object['connected_nodes'], **cell_graph_object['connected_nodes']}
    s_graph, s_positions, s_edge_map = draw_all_line_segments(graph=d_graph,
                                                              positions=d_positions,
                                                              virtual_edge_set=edge_map,
                                                              bounds=outer_bounds,
                                                              already_extended=connected_nodes)
    print(f"\nS edge Map:")
    [print(f"{k} - {v}") for k,v in s_edge_map.items()]
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
        face_vertex_map=cell_graph_object["vertex_map"],
        face_intersection_map=intersection_map)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{output_directory}/graph_6.png")

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
    save_drawn_graph(f"{output_directory}/graph_7.png")

    # Select Subfaces
    print(f"\nSUBFACE SELECTION")
    plane_graph_sub_faces = find_all_subfaces(target_faces=selected_faces,
                                              face_vertex_map=subface_vertex_map,
                                              graph=c_graph)

    # Calculate each subface's centroid
    subface_centroids = get_split_vertex_locations(positions=c_positions,
                                                   target_face_subfaces=plane_graph_sub_faces)

    # Calculate the number of edge crossing induced by connected each subface to all target neighbors
    induced_edge_crossings = calculate_induced_edge_crossings(graph=r_graph,
                                                              positions=r_positions,
                                                              centroids=subface_centroids,
                                                              target_neighbors=target_adjacency)

    # Get Sub_face's edge set
    subfaces_edge_sets = get_face_sub_face_edge_sets(face_sub_cells=plane_graph_sub_faces,
                                                     graph=c_graph)

    print(f"\nsubfaces_edge_sets:")
    print(subfaces_edge_sets)
    print(f"\ninduced edge crossings")
    print(induced_edge_crossings)
    print(f"\nsubfaces:")
    print(plane_graph_sub_faces)
    print(f"\nedge map:")
    print(plane_face_virtual_edge_map)

    # Draw the segment graph
    draw_graph(graph=c_graph, positions=c_positions)
    save_drawn_graph(f"{output_directory}/graph_8.png")

    # Place Copies
    print(f"\nPLACING SPLIT VERTICES")

    # Calculate each subface's centroid
    subface_centroids = get_split_vertex_locations(positions=c_positions,
                                                   target_face_subfaces=plane_graph_sub_faces)

    # Get the number of induced edge crossings in the form of a dictionary of pandas dataframes
    induced_edge_crossing_table = get_edge_crossing_table(induced_edge_crossings=induced_edge_crossings,
                                                          target_neighbors=target_adjacency)

    #
    selected_sub_faces = select_sub_faces(sub_face_tables=induced_edge_crossing_table,
                                          target_faces=selected_faces,
                                          target_vertices=target_adjacency)
    print(f"selected sub_faces: {selected_sub_faces}")

    n_graph, n_positions = place_split_vertices(faces=selected_faces,
                                                selected_sub_faces=selected_sub_faces,
                                                centroids=subface_centroids,
                                                target_vertex=target_vertex,
                                                graph=r_graph,
                                                positions=r_positions)

    # Draw the segment graph
    draw_graph(graph=n_graph, positions=n_positions)
    save_drawn_graph(f"{output_directory}/graph_9.png")

    sys.exit()