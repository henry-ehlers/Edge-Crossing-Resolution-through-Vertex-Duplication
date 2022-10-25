from src.line_segments import *
from src.edge_crossings import *
from src.vertex_splitting import calculate_face_centroid
from src.graph_drawing import *
from src.faces import *

import collections.abc
import networkx as nx
import itertools as it
import pandas as pd
import numpy as np
import copy
import sys


def get_sight_cells_edge_sets(face_sight_cells, graph):
    sight_cell_edge_list = {face: {} for face in face_sight_cells.keys()}
    for face in face_sight_cells.keys():
        sight_cell_edge_list[face].update(get_sight_cell_edges(face_sight_cells[face], graph))
    return sight_cell_edge_list


def get_sight_cell_edges(sight_cells, graph):
    edge_set = {cell: set() for cell in sight_cells}
    for cell in sight_cells:
        cell_edges = get_face_vertex_sequence(cell, graph)
        [edge_set[cell].add(frozenset(edge)) for edge in cell_edges]
    return edge_set


def project_face_sight_lines(edges, vertices, inner_angles, edge_map, graph, positions, bounds, outer):

    # Keep track of the added vertices, and in which edges they were added
    real_nodes = [node for node, real in nx.get_node_attributes(graph, "real").items() if real == 1]
    added_vertices, edge_to_virtual_vertices = [], {}

    # Consider only those vertices whose angle is greater than 180 degrees
    bend_vertices = [key for key in inner_angles.keys() if inner_angles[key] > 180]
    connected_vertices = {vertex: {} for vertex in vertices}
    boundary_edges = nx.get_edge_attributes(graph, "boundary")

    for joint_vertex in bend_vertices:

        if joint_vertex not in real_nodes:
            print("NOT REAL")
            continue

        for connecting_vertex in vertices:
            print(f"\njoint vertex: {joint_vertex} and connecting vertex: {connecting_vertex}")

            # Skip any vertex pair that is a) consists of the same vertex, or b) has already been investigated
            if connecting_vertex == joint_vertex:
                print("THE SAME")
                continue

            if connecting_vertex not in real_nodes:
                print("NOT REAL")
                continue

            # Check whether bend and other vertex can 'see' each other
            is_visible = is_vertex_visible(joint_vertex=joint_vertex,
                                           connecting_vertex=connecting_vertex,
                                           inner_angles=inner_angles,
                                           edge_map=edge_map,
                                           graph=graph,
                                           vertices=vertices,
                                           edges=edges,
                                           positions=positions,
                                           outer=outer)

            # If they cannot see each other, skip to the next pair
            if not is_visible:
                print("NOT VISIBLE")
                continue

            # Extend the sight-line, producing a
            bisected_edge, new_vertex = extend_sight_line(joint_vertex=joint_vertex,
                                                          connecting_vertex=connecting_vertex,
                                                          inner_angles=inner_angles,
                                                          edge_map=edge_map,
                                                          vertices=vertices,
                                                          edges=edges,
                                                          graph=graph,
                                                          positions=positions,
                                                          bounds=bounds,
                                                          outer=outer)

            # Vertices can see one-another, but not produce a legal extension.
            if bisected_edge is None:
                print("NO INTERSECTION")
                continue
            print(f"ADDED VERTEX {new_vertex}")

            # Keep track of what has been added
            added_vertices.append(new_vertex)
            connected_vertices[connecting_vertex][joint_vertex] = (new_vertex, boundary_edges.get(bisected_edge, 0))

            # Add
            if bisected_edge in edge_to_virtual_vertices:
                edge_to_virtual_vertices[bisected_edge].add(new_vertex)
            else:
                edge_to_virtual_vertices[bisected_edge] = {new_vertex}

    return edge_to_virtual_vertices, added_vertices, connected_vertices


def get_faces_sight_cell_incidences(sight_cells, target_vertices, face_edges, face_edge_map, positions):
    # Initialize an empty map of sight cell to incidence
    sight_cell_incidences = {sight_cell: {} for sight_cell in sight_cells}

    # Iterate over all faces in the graph
    for face in sight_cells.keys():
        # Extract all edges in the face, i.e. the virtual edges formed by virtual edge bisection
        face_edge_list = unlist([face_edge_map[edge] for edge in face_edges[face]])

        # Get Incidences of sight cells in current face
        sight_cell_incidences[face].update(get_face_sight_cells_incidences(face_sight_cells=sight_cells[face],
                                                                           face_edge_list=face_edge_list,
                                                                           target_vertices=target_vertices,
                                                                           positions=positions))

    # Return a
    return sight_cell_incidences


def get_outer_face_sight_cell_incidences(sight_cells: {frozenset},
                                         target_vertices: [int],
                                         face_edges: {frozenset: [(int, int)]},
                                         face_edge_map: {frozenset: {frozenset}},
                                         positions):
    [print(f"{key} - {item}") for key, item in face_edge_map.items()]
    # input("+++++++++++++++++++++++++++++++++")

    # Initialize an empty map of sight cell to incidence
    sight_cell_incidences = {sight_cell: set() for sight_cell in sight_cells}

    # Iterate over all faces in the graph
    for cell in sight_cells:
        visibility = [None] * len(face_edges.keys())

        #
        for face_index, face in enumerate(face_edges.keys()):

            # Extract all edges in the face, i.e. the virtual edges formed by virtual edge bisection
            face_edge_list = [tuple(e) for edge in face_edges[face] for e in face_edge_map[frozenset(edge)]]
            print(f"face {face} edge list: {face_edge_list}")
            # Get Incidences of sight cells in current face
            visibility[face_index] = get_sight_cell_incidence(sight_cell_vertices=cell,
                                                              target_vertices=target_vertices,
                                                              real_face_edges=face_edge_list,
                                                              positions=positions)
        sight_cell_incidences[cell] = set.intersection(*visibility)

    # Return a
    return sight_cell_incidences


def get_inner_face_sight_cell_incidences(sight_cells: {frozenset},
                                         target_vertices: [int],
                                         face_edges: {frozenset: [(int, int)]},
                                         face_edge_map: {frozenset: {frozenset}},
                                         positions: {int: np.array}):

    # Initialize an empty map of sight cell to incidence
    sight_cell_incidences = {sight_cell: set() for sight_cell in sight_cells}

    # Iterate over all faces in the graph
    for cell in sight_cells:
        print(f"inner cell: {cell}")

        visibility = [None] * len(face_edges.keys())

        for face_index, face in enumerate(face_edges.keys()):

            # Extract all edges in the face, i.e. the virtual edges formed by virtual edge bisection
            face_edge_list = [tuple(e) for key_edge in face_edges[face] for e in face_edge_map[frozenset(key_edge)]]

            print(f"face edge list: {face_edge_list}")

            # Get Incidences of sight cells in current face
            visibility[face_index] = get_sight_cell_incidence(sight_cell_vertices=cell,
                                                              target_vertices=target_vertices,
                                                              real_face_edges=face_edge_list,
                                                              positions=positions)
        sight_cell_incidences[cell] = set.intersection(*visibility)

    # Return a
    return sight_cell_incidences


def get_face_sight_cells_incidences(face_sight_cells, face_edge_list, target_vertices, positions):
    # If the original face is convex, the sight cell is equivalent
    if len(face_sight_cells) == 1:
        return {face_sight_cells: target_vertices.intersection(get_sorted_face_vertices(face_edge_list))}

    # Iterate over very sight cell in the current face and check visibility to the target vertices
    sight_cell_incidences = {cell: set() for cell in face_sight_cells}
    for sight_cell in face_sight_cells:
        sight_cell_incidences[sight_cell] = get_sight_cell_incidence(sight_cell_vertices=sight_cell,
                                                                     target_vertices=target_vertices,
                                                                     real_face_edges=face_edge_list,
                                                                     positions=positions)
    return sight_cell_incidences


def get_sight_cell_incidence(sight_cell_vertices, target_vertices, real_face_edges, positions):
    """"""

    # Check what targets are already in the current sight cell
    targets_in_cell = sight_cell_vertices.intersection(target_vertices)
    remaining_targets = set(target_vertices) - targets_in_cell

    # Define the sight cell's positions and centroid, as well as initialize incidence set
    sight_cell_positions = [positions[vertex] for vertex in list(sight_cell_vertices)]
    sight_cell_centroid = calculate_face_centroid(sight_cell_positions)
    sight_cell_incidence = set(targets_in_cell)

    # Iterate over all targets that are not already incident to the cell
    for target in remaining_targets:
        sight_line = [sight_cell_centroid, positions[target]]

        # Iterate over all edges and check whether they intersect the line between centroid and target
        is_target_visible = True
        for edge in real_face_edges:
            if target in edge:
                continue

            # Check edge crossings of sight line against edge line extended
            # TODO: floating point issues can cause incorrect results
            edge_line = [positions[vertex] for vertex in edge]
            intersection = line_intersection(sight_line[0], sight_line[1], edge_line[0], edge_line[1])

            # If the two lines intersect, they are not visible
            if intersection is not None:
                is_target_visible = False
                break

        # If no edges intersected the sight line, append it as a visible target
        if is_target_visible:
            sight_cell_incidence.add(target)

    # Return the found set of visible vertices
    return sight_cell_incidence


def merge_face_sight_cells(cells, cells_edge_list, cell_incidences, graph, positions):
    print(f"\nEdge List:")
    [print(f"{key} - {item}") for key, item in cells_edge_list.items()]
    print(f"\nIncidence:")
    [print(f"{key} - {item}") for key, item in cell_incidences.items()]
    for index_a in range(0, len(cells) - 1):
        for index_b in range(index_a + 1, len(cells)):
            cell_a, cell_b = cells[index_a], cells[index_b]

            # Attempt to merge the two cells, return a boolean for success and a (possibly empty) vertex ID
            merge_successful = try_merge_two_sight_cells(cell_a=cell_a,
                                                         cell_b=cell_b,
                                                         cells=cells,
                                                         cells_edge_list=cells_edge_list,
                                                         cell_incidences=cell_incidences,
                                                         graph=graph)

            # If the merge was successful, recurse
            if merge_successful:
                print(f"merged {cell_a} and {cell_b}")

                # draw_graph(graph=graph, positions=positions)
                # save_drawn_graph(f"./{cell_a}_{cell_b}.png")

                # Recurse and repeat merging
                merge_face_sight_cells(cells=cells,
                                       cells_edge_list=cells_edge_list,
                                       cell_incidences=cell_incidences,
                                       graph=graph,
                                       positions=positions)

                # Exit recursion (as the set of cells has changed) and return the removed vertices
                return

    # Final loop, nothing left to merge, return the removed vertices
    return


def try_merge_two_sight_cells(cell_a, cell_b, cells, cells_edge_list, cell_incidences, graph):

    # Determine along which the two cells are to be merged and where the cells are located in the list
    real_edges = [set(edge) for edge, real in nx.get_edge_attributes(G=graph, name="real").items() if real == 1]
    merge_edges = cells_edge_list[cell_a].intersection(cells_edge_list[cell_b])
    incidence_a, incidence_b = cell_incidences[cell_a], cell_incidences[cell_b]
    non_overlapping_incidences = incidence_a ^ incidence_b

    if (non_overlapping_incidences) or (len(merge_edges)) == 0 or (merge_edges in real_edges):
        return False

    # Determine the new cell's vertex set and edge list
    print(f"edge set of {cell_a}: {cells_edge_list[cell_a]}")
    print(f"edge set of {cell_b}: {cells_edge_list[cell_b]}")
    new_edge_set = cells_edge_list[cell_a].union(cells_edge_list[cell_b]) - merge_edges
    new_cell = cell_a.union(cell_b)
    print(f"edge set of NEW {new_cell}: {new_edge_set}")

    # Update the Cell List
    [cells.remove(cell) for cell in [cell_a, cell_b]]
    cells.append(new_cell)

    # Update each Cell's edge list
    [cells_edge_list.pop(cell) for cell in [cell_a, cell_b]]
    cells_edge_list[new_cell] = new_edge_set

    # Update the Cells Incidences
    new_cell_incidence = copy.deepcopy(cell_incidences[cell_a])
    [cell_incidences.pop(cell) for cell in [cell_a, cell_b]]
    cell_incidences[new_cell] = new_cell_incidence

    # Update the graph
    for merge_edge in merge_edges:
        merge_edge = list(merge_edge)
        graph.remove_edge(u=merge_edge[0], v=merge_edge[1])

    # The merge has successfully occurred
    return True


def update_merged_sight_cell_data(face_cells, face_cell_incidences, face_cell_edges, deleted_vertices):
    """

    :param face_cells:
    :param face_cell_incidences:
    :param face_cell_edges:
    :param deleted_vertices: a FROZENSET of vertex ID's
    :return: nothing; the provided three datastructure are modified in place and returned by reference
    """

    # Skip if no vertices were deleted
    if not deleted_vertices:
        return

    # Clear Face's Old set of sight cells
    face_cells.clear()

    # Update Face-Level Information
    old_cells = list(face_cell_incidences.keys())
    for old_cell in old_cells:

        # Check whether the old cell contains any deleted vertices
        if not deleted_vertices.intersection(old_cell):
            continue

        # Create new cell without any of the previously deleted vertices
        new_cell = set(old_cell)
        [new_cell.discard(vertex) for vertex in deleted_vertices]
        new_cell = frozenset(new_cell)

        # Update The List of Cells
        face_cells.add(new_cell)

        # Update the Face Cell Incidences
        face_cell_incidences[new_cell] = face_cell_incidences.pop(old_cell)

        # Update the Dictionary of Cell Edges
        face_cell_edges[new_cell] = face_cell_edges.pop(old_cell)
        edges = copy.copy(face_cell_edges[new_cell])
        [face_cell_edges[new_cell].remove(edge) for edge in edges if deleted_vertices.intersection(edge)]


def match_cell_and_face_incidence(face_incidences, selected_sight_cell_incidences):
    # Indicator whether faces must be reranked (i.e if incidences did not match)
    rerank_faces = False

    # Iterate over all (both selected faces)
    for face in face_incidences.keys():

        # Extract the set of target vertices incident to the current face and its selected cells
        face_incidence, cell_incidences = face_incidences[face], selected_sight_cell_incidences[face]

        # If the number of cells is greater than 1, than no single face contains the target incidence set
        if len(cell_incidences) > 1:
            # Replace the original face with its sight cells
            face_incidences.pop(face)
            face_incidences.update(cell_incidences)

            # Move on to next face and rerank all faces.
            rerank_faces = True
            break

    # Return indicator
    return rerank_faces


def extend_sight_line(joint_vertex, connecting_vertex, inner_angles, edge_map: {frozenset: {frozenset}},
                      vertices, edges, graph, positions, bounds, outer):

    # Calculate intersections of extended line with boundaries in both directions
    bound_intersections = extend_line(positions[joint_vertex], positions[connecting_vertex], bounds)
    distances = [squared_distance(position, positions[connecting_vertex]) for position in bound_intersections]
    closest_intersection_to_joint = bound_intersections[0]
    print(f"\nbound intersections: {bound_intersections}")
    print(f"distances: {distances}")
    print(f"closest: {closest_intersection_to_joint}")
    # input("WHAT?")
    # If vertices are adjacent, they can see one-another; otherwise we must check explicitly
    already_connected = are_vertices_adjacent(joint_vertex, connecting_vertex, graph) \
        or are_vertices_adjacent_virtually(joint_vertex, connecting_vertex, edge_map)

    # If there is only 1 vertex is in the dict of inner angles, skip check
    is_visible = True if len(inner_angles) == 1 \
        else check_vertex_visibility_by_angle(joint_vertex=joint_vertex,
                                              inner_angles=inner_angles,
                                              edges=edges,
                                              vertices=vertices,
                                              positions=positions,
                                              connecting_position=closest_intersection_to_joint,
                                              outer=outer)
    print(f"is visible in extension: {is_visible}")
    # If the hypothetical and observed angle are incompatible, then continue
    if not is_visible:
        return None, None

    # Find the Closest Intersection of the extended line with edges not incident to joint or connecting vertex
    extended_line = (positions[joint_vertex], closest_intersection_to_joint)
    print(f"extended line: {extended_line}")
    # TODO: check also whether edge is part of virtual map => block all those as well
    candidate_edges = [edge for edge in edges if not set(edge).intersection((joint_vertex, connecting_vertex))]
    print(f"candidate edges: {candidate_edges}")
    print(f"positions: {positions}")
    closest_edge, crossing_point = find_closest_edge_intersection(edge_points=extended_line,
                                                                  other_edges=candidate_edges,
                                                                  graph=graph,
                                                                  positions=positions,
                                                                  must_be_real=True)

    # Add Virtual Vertex at Point of Intersection and a virtual edge between it and the origin
    origin_vertex, new_vertex_index = joint_vertex, max(graph.nodes) + 1
    graph.add_node(node_for_adding=new_vertex_index, virtual=1)
    graph.add_edge(u_of_edge=origin_vertex, v_of_edge=new_vertex_index, segment=1)
    positions[new_vertex_index] = np.array(crossing_point)

    # Also add edge between two vertices that were extended (if they do not already have an edge)
    if not already_connected:
        graph.add_edge(u_of_edge=joint_vertex, v_of_edge=connecting_vertex, segment=1)

    # Return a list of added vertices and a map of edges to newly placed virtual vertices
    return closest_edge, new_vertex_index


def is_vertex_visible(joint_vertex, connecting_vertex, inner_angles, edge_map: {frozenset: {frozenset}},
                      graph, vertices, edges, positions, outer):
    # If vertices are neighbors, they can see one-another
    if are_vertices_adjacent(joint_vertex, connecting_vertex, graph) or \
            are_vertices_adjacent_virtually(joint_vertex, connecting_vertex, edge_map):
        return True

    # Check Angle around the Joint Vertex allows for visibility to the connecting vertex
    print(f"inner angles: {inner_angles}")
    angle_visibility = True if len(inner_angles) == 1 else \
        check_vertex_visibility_by_angle(joint_vertex=joint_vertex,
                                         inner_angles=inner_angles,
                                         edges=edges,
                                         vertices=vertices,
                                         positions=positions,
                                         connecting_vertex=connecting_vertex,
                                         outer=outer)
    print(f"angle visibility: {angle_visibility}")

    # If the angle does not allow for visibility, return False
    if not angle_visibility:
        return angle_visibility

    # If the angle (hypothetically) allows for visibility, now check all possible edges that could be in the way
    real_edges = [edge for edge, real in nx.get_edge_attributes(G=graph, name="real").items() if real == 1]
    print(f"real edges: {real_edges}")
    bound_edges = [edge for edge, real in nx.get_edge_attributes(G=graph, name="boundary").items() if real == 1]
    print(f"bound edges: {bound_edges}")
    real_edges = [edge for edge in real_edges if edge not in bound_edges]
    print(f"real edges: {real_edges}")
    candidates = [edge for edge in real_edges if (joint_vertex not in edge) and (connecting_vertex not in edge)]
    print(f"possible crossing edgesS: {candidates}")

    # Check Crossing visibility
    crossing_visibility = check_vertex_visibility_by_crossing(vertex_a=joint_vertex,
                                                              vertex_b=connecting_vertex,
                                                              candidate_edges=candidates,
                                                              positions=positions)
    print(f"crossing visibility: {crossing_visibility}")

    # If the line segment between joint and connecting vertex crosses and edge, return False
    if not crossing_visibility:
        return crossing_visibility

    # All Conditions met for Visibility, so return True
    return True


def check_vertex_visibility_by_angle(joint_vertex, inner_angles, edges, vertices, positions, outer,
                                     connecting_vertex = None, connecting_position = None):

    # Ensure that either a vertex or a position has been provided
    assert (connecting_vertex is not None) or (connecting_position is not None), \
        "Specify a Connecting Point when checking its visibility by angle"

    # Get points for new angle calculation
    joint_index = [index for index in range(0, len(vertices)) if vertices[index] == joint_vertex][0]
    ref_vertex_a, ref_vertex_b = vertices[(joint_index + 1) % len(vertices)], vertices[joint_index - 1]
    debug_angle = calculate_outer_angle(positions[ref_vertex_a], positions[joint_vertex], positions[ref_vertex_b]) \
        if outer else calculate_inner_angle(positions[ref_vertex_a], positions[joint_vertex], positions[ref_vertex_b])

    # Get the Angle of the Joint against which we are comparing the new incoming angle:
    observed_angle = inner_angles[joint_vertex]
    assert observed_angle == debug_angle, \
        f"Angle Calculation is incorrect: Debug ({debug_angle}) and Observed ({observed_angle})"

    # Calculate Hypothetical Angle
    connecting_position = connecting_position if (connecting_position is not None) else positions[connecting_vertex]
    hypothetical_angle = calculate_outer_angle(positions[ref_vertex_a], positions[joint_vertex], connecting_position) \
        if outer else calculate_inner_angle(positions[ref_vertex_a], positions[joint_vertex], connecting_position)

    # If the angle between Vertex A and B is larger than between Vertex A and its neighbors, Vertex B is not visible
    return hypothetical_angle < observed_angle


def check_vertex_visibility_by_crossing(vertex_a, vertex_b, candidate_edges, positions):
    # Get Points of Vertices A and B
    point_a, point_b = positions[vertex_a], positions[vertex_b]

    # For each candidate edge, check crossing between it and (Vertex A, Vertex B)
    for edge in candidate_edges:
        print(f"edge: {edge}")

        # If Lines intersect, Vertices A and B cannot see one-another
        point_c, point_d = positions[edge[0]], positions[edge[1]]
        if line_intersection(point_a, point_b, point_c, point_d) is not None:
            print("INTERSECTION")
            return False

    # If no real edges intersected
    return True


def merge_cells_wrapper(face_sight_cells: {frozenset},
                        cell_incidences: {frozenset: frozenset},
                        cells_edge_map: {frozenset: {frozenset}},
                        ordered_cell_edges: {frozenset: [(int, int)]},
                        positions,
                        graph):
    """"""

    #
    cells_edge_list = {cell: [] for cell in ordered_cell_edges.keys()}
    for cell in ordered_cell_edges.keys():
        cells_edge_list[cell] = set([frozenset(edge) for edge in ordered_cell_edges[cell]])
    print(f"cell_edge_list:")
    print(cells_edge_list)
    # input("CHECK")

    # Try Merging Cells in non-convex face
    face_sight_cells = list(face_sight_cells)  # Reminder: cast is done to enable more efficient indexed looping
    merge_face_sight_cells(cells=face_sight_cells,
                           cells_edge_list=cells_edge_list,
                           cell_incidences=cell_incidences,
                           graph=graph,
                           positions=positions)
    face_sight_cells = set(face_sight_cells)

    # Update the face's cells, their incidents, and edges based on deleted vertices
    cell_vertex_map = update_merged_sight_cells(sight_cells=face_sight_cells,
                                                cell_incidences=cell_incidences,
                                                edge_map=cells_edge_map,
                                                graph=graph)

    print(f"\nPOST MERGE --------------------------------------------------------------------------")
    print(f"\nEdge List:")
    [print(f"{key} - {item}") for key, item in cells_edge_list.items()]
    print(f"\nEdge List:")
    [print(f"{key} - {item}") for key, item in cells_edge_list.items()]
    print(f"\nCells:")
    [print(f"{cell}") for cell in face_sight_cells]
    # input("pass?")
    #
    ordered_cell_edges = {cell: [] for cell in face_sight_cells}
    for cell in face_sight_cells:
        ordered_cell_edges[cell] = get_ordered_edges(edges=[tuple(edge) for edge in list(cells_edge_list[cell])])

    print(f"ordered_cell_edges:")
    print(ordered_cell_edges)
    # input("CHECK")

    # Return updated sight cells, incidences, and edge map
    return face_sight_cells, ordered_cell_edges, cell_vertex_map


def update_merged_sight_cells(sight_cells: {frozenset},
                              cell_incidences: {frozenset: frozenset},
                              edge_map: {frozenset: {frozenset}},
                              graph):

    # Replace all vertices which are virtual, have a degree of 2, and connected only to virtual vertices
    virtual_nodes = nx.get_node_attributes(graph, "virtual")
    cell_vertex_map = {copy.deepcopy(cell): copy.deepcopy(cell) for cell in sight_cells}

    # Iterate over all vertices to determine whether it must be removed
    for node in list(graph.nodes()):

        # Get the current node's edges
        node_edges = list(graph.edges(node))

        # If the node has a degree of more than 2, it must remain
        if len(node_edges) != 0:
            continue

        # If the node is not virtual, do not consider it
        if node not in virtual_nodes:
            continue

        # Remove Edge from edge map and vertex from sight cell
        remove_elements_from_dictionary_frozenset_key(element=node,
                                                      dictionary=cell_incidences)
        remove_elements_from_dictionary_frozenset_key(element=node,
                                                      dictionary=cell_vertex_map)
        # remove_vertex_from_sight_cell(vertex=node,
        #                               sight_cells=sight_cells)

    return cell_vertex_map


def update_sight_cell_graph(sight_cells, cell_incidences, edge_map, graph, positions):

    # Replace all vertices which are virtual, have a degree of 2, and connected only to virtual vertices
    virtual_nodes = nx.get_node_attributes(graph, "virtual")

    # Iterate over all vertices to determine whether it must be removed
    for node in list(graph.nodes()):

        # Get the current node's edges
        node_edges = list(graph.edges(node))

        # If the node has a degree of more than 2, it must remain
        if len(node_edges) > 2 or len(node_edges) == 1:
            continue

        # If the node is not virtual, do not consider it
        if node not in virtual_nodes:
            continue

        # Replace removed vertex's edges with single one
        if len(node_edges) == 2:

            # Create new edge
            vertex_a, vertex_b = set(node_edges[0]) - {node}, set(node_edges[1]) - {node}
            new_edge = (vertex_a.pop(), vertex_b.pop())

            # Replace old edge in map, remove old edge from graph, and add new edge to graph
            replace_edges_with_replacement(edges=node_edges, replacement=new_edge, edge_map=edge_map)
            graph.add_edge(u_of_edge=new_edge[0], v_of_edge=new_edge[1])
            graph.remove_edges_from(node_edges)

        # Remove Vertex from graph and positions list
        graph.remove_node(node)
        positions.pop(node)

        # Remove Edge from edge map and vertex from sight cell
        remove_elements_from_dictionary_frozenset_key(element=node, dictionary=cell_incidences)
        remove_edges_with_vertex(vertex=node, edge_map=edge_map)
        remove_vertex_from_sight_cell(vertex=node, sight_cells=sight_cells)

    # Remove Edges which do not have any virtual maps
    [edge_map.pop(key_edge) for key_edge in list(edge_map.keys()) if (len(edge_map.get(key_edge)) == 0)]

    return sight_cells, cell_incidences, edge_map


def remove_elements_from_dictionary_frozenset_key(element, dictionary):

    #
    for key in list(dictionary.keys()):
        if not key.intersection({element}):
            continue

        #
        new_key = set(key)
        new_key.remove(element)
        dictionary[frozenset(new_key)] = copy.deepcopy(dictionary[key])
        dictionary.pop(key)


def remove_edges_with_vertex(vertex, edge_map):
    for mapped_edge in list(edge_map.keys()):
        if vertex in mapped_edge:
            edge_map.pop(mapped_edge)
            continue
        for sub_edge in list(edge_map[mapped_edge]):
            if vertex not in sub_edge:
                continue
            edge_map[mapped_edge].remove(sub_edge)


def replace_edges_with_replacement(edges, replacement, edge_map):

    # Iterate over all edges in the map
    for mapped_edge in edge_map.keys():

        # Turn all edges into frozensets and lists thereof into sets
        mapped_edge_set = set([frozenset(list(edge)) for edge in edge_map[mapped_edge]])
        remaining_edge_sets = [frozenset(list(edge)) for edge in edges]

        # If both edges (and they should) are in a mapped edge's edge set, replcae it
        if all([edge in mapped_edge_set for edge in remaining_edge_sets]):

            # Remove the found edges from the mapped edge's edge set
            [mapped_edge_set.remove(edge) for edge in remaining_edge_sets]

            # Add the replacement edge
            mapped_edge_set.add(replacement)

            # Replace the mapped edge's edge list with the altered set
            edge_map[mapped_edge] = [tuple(edge) for edge in mapped_edge_set]
            return


def remove_vertex_from_sight_cell(vertex, sight_cells):
    for sight_cell in copy.deepcopy(sight_cells):
        if not sight_cell.intersection({vertex}):
            continue

        sight_cells.remove(sight_cell)
        new_cell = set(sight_cell)
        new_cell.remove(vertex)
        new_cell = frozenset(new_cell)
        sight_cells.add(new_cell)


def update_sight_line_graph(face_vertices: [int],
                            edge_to_virtual_vertices,
                            graph,
                            positions,
                            outer=False):

    print(f"\nface vertices in edge crossing check: {face_vertices}")

    # Remove edges which have been intersected, and replace them with ordered virtual edges
    virtual_edge_map = add_virtual_edges(graph=graph,
                                         positions=positions,
                                         edge_to_virtual_vertex=edge_to_virtual_vertices)  # TODO
    print("virtual edge map:")
    [print(f"{k} - {v}") for k, v in virtual_edge_map.items()]
    remove_edges(graph, edge_to_virtual_vertices.keys())
    print(graph.edges)
    # input("...")

    # Locate Edge Crossings and Faces in Subgraph
    face_positions = {key: positions.get(key) for key in face_vertices}
    print(f"face positions: {face_positions}")
    face_graph = nx.Graph(graph.subgraph(nodes=face_vertices))
    print(f"face graph vertices: {face_graph.nodes()}")
    # input(f"sight line face graph is frozen: {nx.is_frozen(face_graph)}")
    # Find remaining edge crossings (between placed line-segments) and replace them with virtual vertices
    face_edge_crossings, vertex_crossings = locate_edge_crossings(face_graph, face_positions)
    print(f"\nedge crossings: {face_edge_crossings}")
    # input("...")
    if face_edge_crossings:
        virtual_edges = planarize_graph(face_graph, face_positions, face_edge_crossings)
        print(f"virtual edges: {virtual_edges}")
        additional_edges = {frozenset(k): set([frozenset(v) for v in vs]) for k, vs in virtual_edges.items() if vs}
        print(f"non_empty edges: {additional_edges}")
        virtual_edge_map.update(additional_edges)
        graph.update(face_graph)
        positions.update(face_positions)
        for virtual_edge, edges in virtual_edges.items():
            print(f"virtual edge {virtual_edge} - {edges}")
            if len(edges) == 0: continue
            edge_to_be_removed = tuple(virtual_edge)
            graph.remove_edge(u=edge_to_be_removed[0],
                              v=edge_to_be_removed[1])

    # Define Sight Cells, i.e. faces
    print(f"\nUPDATING SIGHT LINE GRAPH:")
    print(graph.edges)
    print()
    print(face_positions)
    sight_cells = find_inner_faces(graph=face_graph, positions=positions) if not outer else None
    if sight_cells is not None:
        print(f"\nplanarizing graph cycles:")
        [print(cell) for cell in sight_cells]
        # input("...")

    # Return Sight Cells
    return sight_cells, virtual_edge_map


def unlist(nested_list):
    return list(it.chain.from_iterable(nested_list))


def project_additional_sight_lines(edges, origin_vertices, origin_angles, target_vertices, target_angles,
                                   edge_map: {frozenset: {frozenset}},
                                   graph, positions, bounds, outer=True):


    real_nodes = [node for node, real in nx.get_node_attributes(graph, "real").items() if real == 1]

    # Keep track of the added vertices, and in which edges they were added
    added_vertices, edge_to_virtual_vertices = [], {}

    # Consider only those vertices whose angle is greater than 180 degrees
    origin_joint_vertices = [key for key in origin_angles.keys() if origin_angles[key] > 180]
    target_joint_vertices = [key for key in target_angles.keys() if target_angles[key] > 180]

    #
    connected_vertices = {vertex: {} for vertex in origin_joint_vertices}
    boundary_edges = nx.get_edge_attributes(graph, "boundary")

    #
    for joint_vertex in origin_joint_vertices:
        print(f"\nJoint Vertex: {joint_vertex}")

        if joint_vertex not in real_nodes:
            continue

        for target_vertex in target_joint_vertices:
            print(f"target vertex: {target_vertex}")

            if target_vertex not in real_nodes:
                continue

            # Check whether bend and other vertex can 'see' each other
            print(f"edges: {edges}")
            is_visible = is_vertex_visible(joint_vertex=joint_vertex,
                                           connecting_vertex=target_vertex,
                                           inner_angles=origin_angles,
                                           edge_map=edge_map,
                                           graph=graph,
                                           vertices=origin_vertices,
                                           edges=edges,
                                           positions=positions,
                                           outer=outer)

            # If they cannot see each other, skip to the next pair
            if not is_visible:
                print(f"IS NOT VISIBLE")
                continue
            print("IS VISIBLE")

            # Extend the sight-line, producing a
            bisected_edge, new_vertex = extend_sight_line(joint_vertex=target_vertex,
                                                          connecting_vertex=joint_vertex,
                                                          inner_angles=target_angles,
                                                          edge_map=edge_map,
                                                          vertices=target_vertices,
                                                          edges=edges,
                                                          graph=graph,
                                                          positions=positions,
                                                          bounds=bounds,
                                                          outer=outer)

            # Vertices can see one-another, but not produce a legal extension.
            if bisected_edge is None:
                continue

            # Keep track of what has been added
            added_vertices.append(new_vertex)

            #
            connected_vertices[joint_vertex][target_vertex] = (new_vertex, boundary_edges.get(bisected_edge, 0))

            # Add
            if bisected_edge in edge_to_virtual_vertices:
                edge_to_virtual_vertices[bisected_edge].add(new_vertex)
            else:
                edge_to_virtual_vertices[bisected_edge] = {new_vertex}

    #
    return edge_to_virtual_vertices, added_vertices, connected_vertices


def add_boundary_to_graph(bounds, graph, positions, offset=0.2, largest_index=None):

    # Specify
    start_index = largest_index + 1 if largest_index is not None else max(graph.nodes) + 1

    # Define the labels of vertices and their edges
    bound_vertices = list(range(start_index, start_index + len(bounds)))
    bound_edges = [(bound_vertices[ind], bound_vertices[(ind + 1) % len(bounds)]) for ind in range(0, len(bounds))]

    # Update Graph, Edges, and Vertex Positions
    for index in range(0, len(bounds)):
        new_position = np.array(bounds[index]) - np.sign(np.array(bounds[index])) * np.array([offset, offset])
        positions[bound_vertices[index]] = new_position
    graph.add_nodes_from(bound_vertices, corner=1)
    graph.add_edges_from(bound_edges, real=1, boundary=1)

    # Return the added vertex labels and edges
    return bound_vertices, bound_edges


def find_inner_face_sight_cells(inner_faces, ordered_face_edges, graph, positions,
                                bounds=((-1, -1), (-1, 1), (1, 1), (1, -1))):

    # Create lists of vertices and edges that define the outer face
    all_face_edges = unlist([ordered_face_edges.get(face) for face in inner_faces if len(face) > 1])
    face_edge_map = {frozenset(edge): {frozenset(edge)} for edge in all_face_edges}
    connected_vertex_map = dict()

    # Iterate over all faces
    selected_face_list = list(inner_faces)
    for face_index, face in enumerate(selected_face_list):

        # Project sight-lines within the currently selected face against itself
        edge_to_virtual_vertices, added_vertices, connected_vertices = project_face_against_self(
            face,
            ordered_face_edges,
            face_edge_map,
            graph,
            positions,
            bounds,
            outer=False)

        if edge_to_virtual_vertices is None:
            continue

        # Update Graph and Virtual Edge Map with New added vertices
        print(f"\nconnected WITHIN FACE {face}: {connected_vertices}")
        update(connected_vertex_map, connected_vertices)
        update_graph_and_virtual_edge_map(face_edge_map,
                                          edge_to_virtual_vertices,
                                          graph,
                                          positions,
                                          outer=False)

        # DEBUG: Draw Initial Embedding
        draw_graph(graph=graph, positions=positions)
        save_drawn_graph(f"./graph_{face}.png")

    # Identify all faces (i.e. sight cells in outer face)
    print(f"\nFIND SIGHT CELLS OF INNER FACES")
    sight_cells, ordered_vertices, ordered_edges = find_inner_faces(graph=graph, positions=positions)
    # print(f"\nINNER SIGHT CELLS")
    # print(f"cells:")
    # [print(cell) for cell in sight_cells]
    # input("...")
    #
    # print(f"ordered cell vertices:")
    # [print(f"cell: {cell} - {vertices}") for cell, vertices in ordered_vertices.items()]
    # input("...")
    #
    # print(f"ordered cell vertices:")
    # [print(f"cell: {cell} - {vertices}") for cell, vertices in ordered_edges.items()]
    # input("...")

    # Return the identified sight cells and the subgraph
    return sight_cells, ordered_vertices, ordered_edges, face_edge_map, connected_vertex_map


def get_clockwise_face_vertices(face, ordered_face_edges, face_edge_map: {frozenset: {frozenset}}, positions,
                                original=False):
    if len(face) == 1:
        return list(face)
    if not original:
        face_edges = unlist([face_edge_map.get(edge) for edge in ordered_face_edges[face]])  # TODO: edges busted
        face_edges = [tuple(edge) for edge in face_edges]
        print(f"mapped face edges of face {face} = {face_edges}")
    else:
        face_edges = [edge for edge in ordered_face_edges[face]]
    face_vertices = get_sorted_face_vertices(face_edges, is_sorted=False)
    if calculate_face_signed_area(face_vertices, positions) < 0:
        face_vertices = list(reversed(face_vertices))
    return face_vertices


def update_graph_and_virtual_edge_map(face_edge_map: {frozenset: {frozenset}}, edge_to_virtual_vertices, graph, positions, outer=False):

    # Update List of Vertices with added vertices and the set of edges which define the face
    face_vertices = graph.nodes

    # Update the Graph and Its Positions
    cells, virtual_edge_map = update_sight_line_graph(face_vertices=face_vertices,
                                                      edge_to_virtual_vertices=edge_to_virtual_vertices,
                                                      graph=graph,
                                                      positions=positions,
                                                      outer=outer)

    #
    print(f"virtual edge map:")
    [print(f"{key} - {item}") for key, item in virtual_edge_map.items()]

    # Update the map of real to virtual edge maps
    deep_update_of_virtual_edge_map(complete_map=face_edge_map,
                                    new_map=virtual_edge_map)

    return cells


def project_face_against_self(face: frozenset,
                              ordered_face_edges: {frozenset: [(int, int)]},
                              face_edge_map: {frozenset: {frozenset}},
                              graph,
                              positions,
                              bounds,
                              outer=False):

    # Get the set of vertices which define the current face
    face_vertices = get_clockwise_face_vertices(face, ordered_face_edges, face_edge_map, positions, original=True)
    print(f"face vertices: {face_vertices}")

    # Calculate Inner Angles to identify joint vertices
    if outer:
        face_angles = calculate_face_outer_angles(counter_clockwise_face_vertices=face_vertices,
                                                  positions=positions)
    else:
        face_angles = calculate_face_inner_angles(counter_clockwise_face_vertices=face_vertices,
                                                  positions=positions)
        if is_convex(face_angles):
            return None, None, None

    print(f"face angles in project_face_against_self: {face_angles}")

    # Replace original edges with their virtual counterparts
    candidate_edges = [tuple(edge) for edge_key in face_edge_map.keys() for edge in face_edge_map.get(edge_key)]
    print(f"candidate edges: {candidate_edges}")

    # Extend all sight lines, add new vertices where necessary, and keep track of edge bisection
    edge_to_virtual_vertices, added_vertices, connected_vertices = project_face_sight_lines(
        edges=candidate_edges,
        vertices=face_vertices,
        inner_angles=face_angles,
        edge_map=face_edge_map,
        graph=graph,
        positions=positions,
        bounds=bounds,
        outer=outer)

    [connected_vertices.pop(v) for v in connected_vertices.keys() if connected_vertices.get(v, None) is None]

    # Return
    return edge_to_virtual_vertices, added_vertices, connected_vertices


def project_outer_face_against_singleton(face: frozenset,
                                         other_face: frozenset,
                                         edge_map: {frozenset: {frozenset}},
                                         ordered_face_edges: {frozenset: [(int, int)]},
                                         positions,
                                         graph,
                                         bounds):
    assert(len(other_face) == 1), \
        f"ERROR: Passed face {other_face} of length {len(other_face)} into singleton projection!"

    # Replace original edges with their virtual counterparts
    candidate_edges = [tuple(edge) for edge_key in edge_map.keys() for edge in edge_map.get(edge_key)]
    other_face_vertices = list(other_face)
    other_face_angles = {vertex: 360 for vertex in other_face_vertices}

    # Get Vertices and ensure they are listed in counter-clockwise order
    face_vertices = get_clockwise_face_vertices(face, ordered_face_edges, edge_map, positions, original=True)
    face_angles = calculate_face_outer_angles(face_vertices, positions)

    #
    edge_to_virtual_vertices, added_vertices, connected_vertices = project_additional_sight_lines(
        edges=candidate_edges,
        origin_vertices=face_vertices,
        origin_angles=face_angles,
        target_vertices=other_face_vertices,
        edge_map=edge_map,
        target_angles=other_face_angles,
        graph=graph,
        positions=positions,
        bounds=bounds,
        outer=True)

    #
    [connected_vertices.pop(v) for v in connected_vertices.keys() if connected_vertices.get(v, None) is None]

    return edge_to_virtual_vertices, added_vertices, connected_vertices


def project_outer_face_against_non_face():
    pass


def project_outer_face_against_another_face(face, other_face, edge_map: {frozenset: {frozenset}}, ordered_face_edges, positions, graph, bounds):

    # Replace original edges with their virtual counterparts
    candidate_edges = [tuple(edge) for edge_key in edge_map.keys() for edge in edge_map.get(edge_key)]
    other_face_vertices = get_clockwise_face_vertices(
        other_face, ordered_face_edges, edge_map, positions, original=True)
    other_face_angles = calculate_face_outer_angles(other_face_vertices, positions)

    # Get Vertices and ensure they are listed in counter-clockwise order
    face_vertices = get_clockwise_face_vertices(face, ordered_face_edges, edge_map, positions, original=True)
    face_angles = {vertex: 360 for vertex in face_vertices} if len(face) == 1 else \
        calculate_face_outer_angles(face_vertices, positions)

    #
    edge_to_virtual_vertices, added_vertices, connected_vertices = project_additional_sight_lines(
        edges=candidate_edges,
        origin_vertices=face_vertices,
        origin_angles=face_angles,
        target_vertices=other_face_vertices,
        edge_map=edge_map,
        target_angles=other_face_angles,
        graph=graph,
        positions=positions,
        bounds=bounds,
        outer=True)

    #
    [connected_vertices.pop(v) for v in connected_vertices.keys() if connected_vertices.get(v, None) is None]

    #
    return edge_to_virtual_vertices, added_vertices, connected_vertices


def get_subgraph(nodes, edges, graph, positions):
    sub_graph, sub_positions = copy.deepcopy(graph), copy.deepcopy(positions)
    sub_graph.remove_nodes_from([n for n in graph if n not in set(nodes)])
    [sub_positions.pop(n, None) for n in graph if n not in nodes]
    sub_graph.remove_edges_from([e for e in sub_graph.edges if set(e) not in [set(edge) for edge in edges]])
    return sub_graph, sub_positions


def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def find_outer_face_sight_cells(selected_faces: {frozenset},
                                ordered_face_edges: {frozenset: [(int, int)]},
                                graph,
                                positions: {int: np.array},
                                is_cycle: [bool],
                                bounds):
    # TODO: a bisected REAL edge will not be extended since we are looking up the original edge sets, which don't

    # Create lists of vertices and edges that define the outer face
    all_face_edges = unlist([ordered_face_edges.get(face) for face in selected_faces if len(face) > 1])
    outer_face_vertices = list(frozenset().union(*selected_faces))

    print(f"\nselected faces: {selected_faces}")
    print(f"\nouter face vertices: {outer_face_vertices}")

    # Create A graph consisting only of outer-face defining vertices and edges
    outer_graph, outer_positions = get_subgraph(outer_face_vertices, all_face_edges, graph, positions)

    # Define the embedding boundary and the face edge map
    bound_vertices, bound_edges = add_boundary_to_graph(bounds=bounds,
                                                        graph=outer_graph,
                                                        positions=outer_positions,
                                                        largest_index=max(graph.nodes))
    face_edge_map = {frozenset(edge): {frozenset(edge)} for edge in all_face_edges + bound_edges}
    connected_vertex_map = dict()

    # Iterate over all faces
    selected_face_list = list(selected_faces)
    for face_index, face in enumerate(selected_face_list):
        if is_cycle[face]:

            # Add additional candidate edges if we are dealing with the outer face
            other_faces = copy.copy(selected_faces)
            other_faces.remove(face)

            # Project sight-lines within the currently selected face against itself
            edge_to_virtual_vertices, added_vertices, connected_vertices = project_face_against_self(
                face=face,
                ordered_face_edges=ordered_face_edges,
                face_edge_map=face_edge_map,
                graph=outer_graph,
                positions=outer_positions,
                bounds=bounds,
                outer=True)

            # Update Graph and Virtual Edge Map with New added vertices
            print(f"\nconnected WITHIN FACE: {connected_vertices}")
            update(connected_vertex_map, connected_vertices)
            print(f"map: {connected_vertex_map}")
            update_graph_and_virtual_edge_map(face_edge_map=face_edge_map,
                                              edge_to_virtual_vertices=edge_to_virtual_vertices,
                                              graph=outer_graph,
                                              positions=outer_positions,
                                              outer=True)

            # DEBUG: Draw Initial Embedding
            draw_graph(graph=outer_graph, positions=outer_positions)
            save_drawn_graph(f"./graph_{face}.png")

        # Iterate over all other faces that form the outer face
        for other_face in selected_faces:
            if other_face == face:
                continue

            # Check if the Other face is a cycle or not -> check angle visibilities + edge crossing visibilities
            if is_cycle[other_face]:
                edge_to_virtual_vertices, added_vertices, connected_vertices = project_outer_face_against_another_face(
                    face,
                    other_face,
                    face_edge_map,
                    ordered_face_edges,
                    outer_positions,
                    outer_graph,
                    bounds)
            elif not is_cycle[other_face] and len(other_face) > 1:
                # TODO: if other face is not a cycle -> special visiblility checking of ONLY edge intersections
                continue
            elif not is_cycle[other_face] and len(other_face) == 1:
                edge_to_virtual_vertices, added_vertices, connected_vertices = project_outer_face_against_singleton(
                    face,
                    other_face,
                    face_edge_map,
                    ordered_face_edges,
                    outer_positions,
                    outer_graph,
                    bounds)
            else:
                sys.exit("Non Accounted for Constellation of face!")

            # Update Graph and Virtual Edge Map with New added vertices
            print(f"\nconnected OUTSIDE FACE: {connected_vertices}")
            update(connected_vertex_map, connected_vertices)
            print(f"map: {connected_vertex_map}")
            update_graph_and_virtual_edge_map(face_edge_map=face_edge_map,
                                              edge_to_virtual_vertices=edge_to_virtual_vertices,
                                              graph=outer_graph,
                                              positions=outer_positions,
                                              outer=True)

            # Draw Initial Embedding
            draw_graph(graph=outer_graph, positions=outer_positions)
            save_drawn_graph(f"./graph_{face}+{other_face}.png")

    print(f"\nfinal map: {connected_vertex_map}")
    draw_graph(graph=outer_graph,
               positions=outer_positions)
    save_drawn_graph(f"./final_outer.png")

    # Identify all faces (i.e. sight cells in outer face)
    sight_cells, ordered_cell_nodes, ordered_cell_edges = find_inner_faces(graph=outer_graph,
                                                                           positions=outer_positions)

    # Remove the original set of faces that defined the outer face
    print(sight_cells)
    [sight_cells.remove(cell) for cell in copy.copy(sight_cells) for face in selected_faces
     if face.issubset(cell) and len(face) >= 2]

    # Return the identified sight cells and the subgraph
    return sight_cells, ordered_cell_edges, face_edge_map, connected_vertex_map, outer_graph, outer_positions


def merge_edge_map(old_edge_map: {frozenset: {frozenset}}, new_edge_map: {frozenset: {frozenset}}):
    old_edge_keys, new_edge_keys = set(list(old_edge_map.keys())), set(list(new_edge_map.keys()))
    complete_edge_keys = old_edge_keys.union(new_edge_keys)
    complete_edge_map = {edge: {edge} for edge in complete_edge_keys}
    deep_update_of_virtual_edge_map(complete_edge_map, old_edge_map)
    deep_update_of_virtual_edge_map(complete_edge_map, new_edge_map)
    return complete_edge_map


def merge_connected_nodes(list_connected_nodes):
    all_connected_nodes = {}
    for connected_nodes in list_connected_nodes:
        for connected_vertex in connected_nodes.keys():
            if connected_vertex not in all_connected_nodes.keys():
                all_connected_nodes[connected_vertex] = {}
            all_connected_nodes[connected_vertex].update(connected_nodes[connected_vertex])
    vertices = list(all_connected_nodes.keys())
    [all_connected_nodes.pop(v) for v in vertices if len(all_connected_nodes.get(v)) == 0]
    return all_connected_nodes


def deep_update_of_virtual_edge_map(complete_map: {frozenset: {frozenset}}, new_map: {frozenset: {frozenset}}):

    print("-------------------------------------------------------------")
    [print(f"{key} - {item}") for key, item in complete_map.items()]
    print("-------------------------------------------------------------")
    [print(f"{key} - {item}") for key, item in new_map.items()]
    # input("CHECK STATUS")
    for intersected_edge, virtual_edges in new_map.items():
        mapped = [real_edge for real_edge in complete_map.keys() if intersected_edge in complete_map[real_edge]]

        # Ensure that the new virtual edge mapped only to a single existing edge
        # assert(len(mapped) <= 1), \
        #     f"Virtual Edge Map found too many mappings of new virtual edge {intersected_edge} to {mapped}"

        # Replace already virtual edge with a new list of virtual edges
        if len(mapped) >= 1:

            for existing_edge in mapped:
                # Extract the intersected edge and identify where it is in the complete map
                # existing_edge = mapped.pop(0)  # can only be of length 1
                print(f"mapped to existing edge: {existing_edge} @ {intersected_edge}")
                complete_map[existing_edge].remove(intersected_edge)  # remove virtual edge from list of complete map\
                print(f"virtual edges to be added: {virtual_edges}")
                print(f"original edges in the set: {complete_map[existing_edge]}")
                complete_map[existing_edge] = complete_map[existing_edge].union(virtual_edges)
                print(f"edges after merger: {complete_map[existing_edge]}")

        # Connection between two vertices which were previously unconnected
        else:

            # Update the complete edge map with a new dictionary entry and list of virtual edges
            print(f"adding {intersected_edge} -> {virtual_edges}")
            complete_map.update({intersected_edge: virtual_edges})


def is_convex(inner_angles):
    return all(angle <= 180.0 for angle in inner_angles.values())


def are_vertices_adjacent(vertex_a, vertex_b, graph):
    # TODO: some edges end up missing a 'real' tag somehow
    adjacent = True
    try:
        graph.edges[vertex_a, vertex_b]
    except KeyError:
        adjacent = False
    if not adjacent:
        return adjacent
    fields = graph.get_edge_data(vertex_a, vertex_b, default=0)
    if not fields:
        return True
    if not (fields.get("real", 0) == 1):
        return False
    return adjacent


def are_vertices_adjacent_virtually(vertex_a, vertex_b, edge_map: {frozenset: {frozenset}}):  # TODO: rename function to are_real_vertices...
    print(f"keys; {list(edge_map.keys())}")
    print(f"edge: {frozenset((vertex_a, vertex_b))}")
    return frozenset((vertex_a, vertex_b)) in list(edge_map.keys())
