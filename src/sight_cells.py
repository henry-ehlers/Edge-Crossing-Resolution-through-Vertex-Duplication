from src.line_segments import *
from src.edge_crossings import *
from src.vertex_splitting import calculate_face_centroid
from src.graph_drawing import *
from src.faces import *

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
    added_vertices, edge_to_virtual_vertices = [], {}
    # Consider only those vertices whose angle is greater than 180 degrees
    bend_vertices = [key for key in inner_angles.keys() if inner_angles[key] > 180]

    for joint_vertex in bend_vertices:
        for connecting_vertex in vertices:

            # Skip any vertex pair that is a) consists of the same vertex, or b) has already been investigated
            if connecting_vertex == joint_vertex:
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
                continue

            # Keep track of what has been added
            added_vertices.append(new_vertex)

            # Add
            if bisected_edge in edge_to_virtual_vertices:
                edge_to_virtual_vertices[bisected_edge].add(new_vertex)
            else:
                edge_to_virtual_vertices[bisected_edge] = {new_vertex}

    return edge_to_virtual_vertices, added_vertices


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


def get_outer_face_sight_cell_incidences(sight_cells, target_vertices, face_edges, face_edge_map, positions):
    # Initialize an empty map of sight cell to incidence
    sight_cell_incidences = {sight_cell: set() for sight_cell in sight_cells}

    # Iterate over all faces in the graph
    for cell in sight_cells:
        visibility = [None] * len(face_edges.keys())

        for face_index, face in enumerate(face_edges.keys()):
            # Extract all edges in the face, i.e. the virtual edges formed by virtual edge bisection
            face_edge_list = unlist([face_edge_map[edge] for edge in face_edges[face]])

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


def merge_face_sight_cells(cells, cells_edge_list, cell_incidences, removed_vertices, graph):
    # todo: outer face also contains all faces INSIDE of the graph

    for cell_index_a in range(0, len(cells)):
        for cell_index_b in range(cell_index_a + 1, len(cells)):
            cell_a, cell_b = cells[cell_index_a], cells[cell_index_b]

            # Attempt to merge the two cells, return a boolean for success and a (possibly empty) vertex ID
            merge_successful, removed = try_merge_two_sight_cells(cell_a=cell_a,
                                                                  cell_b=cell_b,
                                                                  cells=cells,
                                                                  cells_edge_list=cells_edge_list,
                                                                  cell_incidences=cell_incidences,
                                                                  graph=graph)

            # If the merge was successful, recurse
            if merge_successful:

                # If a vertex was removed, keep track if its removal
                if removed:
                    [removed_vertices.append(removed_vertex) for removed_vertex in removed]

                # Recurse and repeat merging
                removed_vertices = merge_face_sight_cells(cells=cells,
                                                          cells_edge_list=cells_edge_list,
                                                          cell_incidences=cell_incidences,
                                                          removed_vertices=removed_vertices,
                                                          graph=graph)

                # Exit recursion (as the set of cells has changed) and return the removed vertices
                return removed_vertices

    # Final loop, nothing left to merge, return the removed vertices
    return removed_vertices


def try_merge_two_sight_cells(cell_a, cell_b, cells, cells_edge_list, cell_incidences, graph):

    # Determine along which the two cells are to be merged and where the cells are located in the list
    merge_edges = cells_edge_list[cell_a].intersection(cells_edge_list[cell_b])
    incidence_a, incidence_b = cell_incidences[cell_a], cell_incidences[cell_b]
    non_overlapping_incidences = incidence_a ^ incidence_b

    if non_overlapping_incidences or len(merge_edges) == 0:
        return False, None

    # Determine the new cell's vertex set and edge list
    new_edge_set = cells_edge_list[cell_a].union(cells_edge_list[cell_b] - merge_edges)
    new_cell = cell_a.union(cell_b)

    # Update the Cell List
    cells.append(new_cell)
    [cells.remove(cell_key) for cell_key in [cell_a, cell_b]]

    # Update each Cell's edge list
    cells_edge_list[new_cell] = new_edge_set
    [cells_edge_list.pop(cell_key, None) for cell_key in [cell_a, cell_b]]

    # Update the Cells Incidences
    cell_incidences[new_cell] = copy.deepcopy(cell_incidences[cell_a])
    [cell_incidences.pop(cell_key, None) for cell_key in [cell_a, cell_b]]

    # Update the graph
    common_vertices = []
    for merge_edge in merge_edges:
        update_merge_sight_cell_graph(merge_edge, graph)
        # todo: make common vertex merging work
        # merge_edge = list(merge_edge)
        # common_vertex = [vertex for vertex in merge_edge[0] if vertex in merge_edge[1] and vertex in merge_edge[0]]
        # print(f"common vertex: {common_vertex}")
        # common_vertices.append(common_vertex)

    # The merge has successfully occurred
    return True, common_vertices


def update_merge_sight_cell_graph(merge_edge, graph):
    merge_edge = list(merge_edge)
    graph.remove_edge(u=merge_edge[0], v=merge_edge[1])

    # Identify a vertex connected to deleted edge which is couched between virtual edges. Remove vertex and replace edge
    vertex_edges, common_vertex = unlist([list(graph.edges(v)) for v in merge_edge if len(graph.edges(v)) == 2]), None

    # todo: maybe clean up the graph of non-connected and virtual edges?
    # if vertex_edges:
    #     # Define a new edge which skips the now singleton vertex
    #     new_edge = [vertex for vertex in unlist(vertex_edges) if vertex not in merge_edge]
    #     real_edge = all([get_graph_entity_data(graph.nodes, vertex, "virtual", 0) for vertex in new_edge])
    #     graph.add_edge(u_of_edge=new_edge[0], v_of_edge=new_edge[1], virtual=0 if real_edge else 1)
    #     print(f"Removed Edges {vertex_edges} with {new_edge}")
    #     # Remove the singleton virtual vertex from the graph and cells
    #     common_vertex = set(vertex_edges[0]).intersection(set(vertex_edges[1])).pop()
    #     graph.remove_node(common_vertex)
    #     print(f"REMOVING VERTEX {common_vertex}")
    #     return common_vertex

    return None


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


def find_minimal_sight_cell_set(cell_incidences):
    """

    :param cell_incidences: SUBSET OF ALL CELL INCIDENCES FOR A PARTICULAR FACE
    :return:
    """
    # TODO: implement Anais' thing here instead. we don;t need one task per one worker, but multiple per worker
    # todo: get target vertex list from list of unique vertices in incidences
    # Initiate cost matrix of ones
    # cost_matrix = np.ones(shape=(len(face_cells_incidence), len(target_vertices)), dtype=int)
    #
    # # Store sight cells in list to avoid ordering problems
    # row_names = list(face_cells_incidence.keys())
    # # Iterate over all sight cells and extract their incidences to build cost matrix
    # for sight_cell in row_names:
    #     row_index = row_names.index(sight_cell)
    #     for visible_vertex in face_cells_incidence[sight_cell]:
    #         col_index = target_vertices.index(visible_vertex)
    #         cost_matrix[row_index, col_index] -= 1

    # Find minimal assignment cost
    incidence_number = {cell: len(cell_incidences[cell]) for cell in cell_incidences.keys()}
    incidence_number = dict(sorted(incidence_number.items(), key=lambda item: item[1], reverse=True))
    selected_cells = list(incidence_number)[0:2]  # todo: selection could be arbitrarily long
    return {cell: cell_incidences[cell] for cell in selected_cells}


def get_sorted_sight_cell_incidence_table(sight_cell_incidences):

    # Collect list of cells and their incidences
    incidence_table = []
    [incidence_table.append((cell, sight_cell_incidences.get(cell, None))) for cell in sight_cell_incidences.keys()]

    # Transform list into pandas dataframe
    incidence_table = pd.DataFrame(data=incidence_table, columns=['sight_cell', 'incidence'])
    incidence_table["n_incidence"] = [len(incidence) for incidence in incidence_table["incidence"]]

    # Sort incidence table by the number of vertices incident to each cell
    incidence_table.sort_values(by='n_incidence', ascending=False, inplace=True, ignore_index=True)

    # Return the table
    return incidence_table


def select_sight_cells(cell_incidences, target_vertices):
    """"""

    # Get the incidence table
    incidence_table = get_sorted_sight_cell_incidence_table(cell_incidences)
    print(incidence_table)

    # Iterate over all cells and their incidences
    selected_cells, remaining_targets = {}, set(target_vertices)
    for index, row in incidence_table.iterrows():

        # If no targets remain, return the selection
        if len(remaining_targets) == 0:
            return selected_cells

        # If the current cell's incidence is not yet selected, select it and update remaining targets
        if remaining_targets.intersection(row['incidence']):
            remaining_targets -= row['incidence']
            selected_cells.update({row['sight_cell']:  row['incidence']})

    # TODO: we should never return here; maybe cover this case somehow?
    return selected_cells


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


def extend_sight_line(joint_vertex, connecting_vertex, inner_angles, edge_map, vertices, edges,
                      graph, positions, bounds, outer):

    # Calculate intersections of extended line with boundaries in both directions
    bound_intersections = extend_line(positions[joint_vertex], positions[connecting_vertex], bounds)
    closest_intersection_to_joint = bound_intersections[0]

    # If vertices are adjacent, they can see one-another; otherwise we must check explicitly
    already_connected = are_vertices_adjacent(joint_vertex, connecting_vertex, graph) \
        or are_vertices_adjacent_virtually(joint_vertex, connecting_vertex, edge_map)
    is_visible = check_vertex_visibility_by_angle(joint_vertex=joint_vertex,
                                                  inner_angles=inner_angles,
                                                  edges=edges,
                                                  vertices=vertices,
                                                  positions=positions,
                                                  connecting_position=closest_intersection_to_joint,
                                                  outer=outer)

    # If the hypothetical and observed angle are incompatible, then continue
    if not is_visible:
        return None, None

    # Find the Closest Intersection of the extended line with edges not incident to joint or connecting vertex
    extended_line = (positions[joint_vertex], closest_intersection_to_joint)
    # TODO: check also whether edge is part of virtual map => block all those as well
    candidate_edges = [edge for edge in edges if not set(edge).intersection((joint_vertex, connecting_vertex))]
    closest_edge, crossing_point = find_closest_edge_intersection(
        extended_line, candidate_edges, graph, positions, must_be_real=True)

    # Add Virtual Vertex at Point of Intersection and a virtual edge between it and the origin
    origin_vertex, new_vertex_index = joint_vertex, max(graph.nodes) + 1

    graph.add_node(node_for_adding=new_vertex_index, virtual=1)
    graph.add_edge(u_of_edge=origin_vertex, v_of_edge=new_vertex_index, segment=1, real=0)
    positions[new_vertex_index] = np.array(crossing_point)

    # Also add edge between two vertices that were extended (if they do not already have an edge)
    if not already_connected:
        graph.add_edge(u_of_edge=joint_vertex, v_of_edge=connecting_vertex, segment=1, real=0)

    # Return a list of added vertices and a map of edges to newly placed virtual vertices
    return closest_edge, new_vertex_index


def is_vertex_visible(joint_vertex, connecting_vertex, inner_angles, edge_map, graph, vertices, edges, positions,
                      outer):
    # If vertices are neighbors, they can see one-another
    if are_vertices_adjacent(joint_vertex, connecting_vertex, graph) or \
            are_vertices_adjacent_virtually(joint_vertex, connecting_vertex, edge_map):
        return True

    # Check Angle around the Joint Vertex allows for visibility to the connecting vertex
    angle_visibility = check_vertex_visibility_by_angle(joint_vertex=joint_vertex,
                                                        inner_angles=inner_angles,
                                                        edges=edges,
                                                        vertices=vertices,
                                                        positions=positions,
                                                        connecting_vertex=connecting_vertex,
                                                        outer=outer)

    # If the angle does not allow for visibility, return False
    if not angle_visibility:
        return angle_visibility

    # If the angle (hypothetically) allows for visibility, now check all possible edges that could be in the way
    possible_crossing_edges = [edge for edge in edges if (joint_vertex not in edge) and (connecting_vertex not in edge)]
    crossing_visibility = check_vertex_visibility_by_crossing(vertex_a=joint_vertex,
                                                              vertex_b=connecting_vertex,
                                                              candidate_edges=possible_crossing_edges,
                                                              positions=positions)

    # If the line segment between joint and connecting vertex crosses and edge, return False
    if not crossing_visibility:
        return crossing_visibility

    # All Conditions met for Visibility, so return True
    return True


def check_vertex_visibility_by_angle(joint_vertex, inner_angles, edges, vertices, positions, outer,
                                     connecting_vertex=None, connecting_position=None):

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

        # If Lines intersect, Vertices A and B cannot see one-another
        point_c, point_d = positions[edge[0]], positions[edge[1]]
        if line_intersection(point_a, point_b, point_c, point_d) is not None:
            return False

    # If no real edges intersected
    return True


def merge_cells_wrapper(face_sight_cells, cells_edge_list, cell_incidences, graph):

    # Try Merging Cells in non-convex face
    face_sight_cells = list(face_sight_cells) # todo: this is does not work
    removed_vertices = merge_face_sight_cells(cells=face_sight_cells,
                                              cells_edge_list=cells_edge_list,
                                              cell_incidences=cell_incidences,
                                              removed_vertices=[],
                                              graph=graph)
    print(face_sight_cells)
    # Update the face's cells, their incidents, and edges based on deleted vertices
    # TODO: replace with update_sight_cell_graph()
    # update_merged_sight_cell_data(face_cells=face_sight_cells,
    #                               face_cell_incidences=cell_incidences,
    #                               face_cell_edges=cells_edge_list,
    #                               deleted_vertices=frozenset(removed_vertices))
    return set(face_sight_cells)


def update_sight_cell_graph(sight_cells, edge_map, graph, positions):
    print(f"\nsight cells:")
    [print(cell) for cell in sight_cells]

    # Replace all vertices which are virtual, have a degree of 2, and connected only to virtual vertices
    virtual_nodes = nx.get_node_attributes(graph, "virtual")

    # Iterate over all vertices to determine whether it must be removed
    for node in list(graph.nodes()):

        # Get the current node's edges
        node_edges = list(graph.edges(node))
        print(f"\nnode: {node} with edges {node_edges}")

        # If the vertex is disconnected, remove it
        if len(node_edges) == 0:
            print("SINGLETON")
            remove_vertex_from_sight_cell(vertex=node, sight_cells=sight_cells)
            remove_edges_with_vertex(vertex=node, edge_map=edge_map)
            graph.remove_node(node)
            continue

        # If the node has a degree of more than 2, it must remain
        if len(node_edges) > 2 or len(node_edges) == 1:
            continue

        # If the node is not virtual, do not consider it
        if node not in virtual_nodes:
            continue

        # Get vertices adjacent to current node and connect them with a new edge
        vertex_a, vertex_b = set(node_edges[0]) - {node}, set(node_edges[1]) - {node}
        new_edge = (vertex_a.pop(), vertex_b.pop())
        graph.add_edge(u_of_edge=new_edge[0], v_of_edge=new_edge[1])

        # Remove the original vertex, its position, and its edges
        graph.remove_edges_from(node_edges)
        graph.remove_node(node)
        positions.pop(node)

        # Remove Edge from edge map and vertex from sight cell
        replace_edges_with_replacement(edges=node_edges, replacement=new_edge, edge_map=edge_map)
        remove_vertex_from_sight_cell(vertex=node, sight_cells=sight_cells)

    print(f"\nsight cells:")
    [print(cell) for cell in sight_cells]
    print(f"\nedge map:")
    [print(f"{edge} - {edge_map[edge]}") for edge in edge_map.keys()]  # TODO: THERE ARE REMAINING EDGES DANGLING

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
            print(f"replaced at {mapped_edge} -> {edge_map[mapped_edge]}")
            return


def remove_vertex_from_sight_cell(vertex, sight_cells):
    for sight_cell in copy.deepcopy(sight_cells):
        if not sight_cell.intersection({vertex}):
            continue

        sight_cells.remove(sight_cell)
        new_cell = set(sight_cell)
        new_cell.remove(vertex)
        new_cell = frozenset(new_cell)
        print(f"found {vertex} in {sight_cell} and replaced with {new_cell}")
        sight_cells.add(new_cell)


def merge_all_face_cells(face_sight_cells, face_cell_edge_map, cell_incidences, graph):
    # Iterate over all faces, and attempt to merge their faces
    for face in face_sight_cells.keys():

        # Skip convex faces
        if len(face_sight_cells[face]) == 1:
            continue

        merge_cells_wrapper(face_sight_cells=face_sight_cells.get(face, None),
                            cells_edge_list=face_cell_edge_map.get(face, None),
                            cell_incidences=cell_incidences.get(face, None),
                            graph=graph)


def update_sight_line_graph(edges, face_vertices, edge_to_virtual_vertices, graph, positions, outer=False):

    # Remove edges which have been intersected, and replace them with ordered virtual edges
    virtual_edge_map = add_virtual_edges(graph, positions, edge_to_virtual_vertices)
    remove_edges(graph, edge_to_virtual_vertices.keys())

    # Locate Edge Crossings and Faces in Subgraph
    face_positions = {key: positions.get(key) for key in face_vertices}
    face_graph = nx.Graph(graph.subgraph(nodes=face_vertices))

    # Find remaining edge crossings (between placed line-segments) and replace them with virtual vertices
    face_edge_crossings, vertex_crossings = locate_edge_crossings(face_graph, face_positions)

    if face_edge_crossings:
        face_graph, face_positions, virtual_edges = planarize_graph(face_graph, face_positions, face_edge_crossings)
        non_empty_virtual_edges = {k: v for k, v in virtual_edges.items() if v}
        virtual_edge_map.update(non_empty_virtual_edges)
        graph.update(face_graph)
        positions.update(face_positions)
        [graph.remove_edge(u=edge[0], v=edge[1]) for edge in virtual_edges.keys() if virtual_edges[edge]]

    # Define Sight Cells, i.e. faces
    sight_cells = find_inner_faces(face_graph, positions=positions) if not outer else None

    # Return Sight Cells
    return sight_cells, virtual_edge_map


def unlist(nested_list):
    return list(it.chain.from_iterable(nested_list))


def project_additional_sight_lines(edges, origin_vertices, origin_angles, target_vertices, target_angles, edge_map,
                                   graph, positions, bounds, outer=True):

    # Keep track of the added vertices, and in which edges they were added
    added_vertices, edge_to_virtual_vertices = [], {}

    # Consider only those vertices whose angle is greater than 180 degrees
    origin_joint_vertices = [key for key in origin_angles.keys() if origin_angles[key] > 180]
    target_joint_vertices = [key for key in target_angles.keys() if target_angles[key] > 180]

    for joint_vertex in origin_joint_vertices:
        for target_vertex in target_joint_vertices:

            # Check whether bend and other vertex can 'see' each other
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
                continue

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

            # Add
            if bisected_edge in edge_to_virtual_vertices:
                edge_to_virtual_vertices[bisected_edge].add(new_vertex)
            else:
                edge_to_virtual_vertices[bisected_edge] = {new_vertex}

    return edge_to_virtual_vertices, added_vertices


def add_boundary_to_graph(bounds, graph, positions, offset=0.2):
    # Define the labels of vertices and their edges
    bound_vertices = list(range(max(graph.nodes) + 1, max(graph.nodes) + 1 + len(bounds)))
    bound_edges = [(bound_vertices[ind], bound_vertices[(ind + 1) % len(bounds)]) for ind in range(0, len(bounds))]

    # Update Graph, Edges, and Vertex Positions
    for index in range(0, len(bounds)):
        new_position = np.array(bounds[index]) - np.sign(np.array(bounds[index])) * np.array([offset, offset])
        positions[bound_vertices[index]] = new_position
    graph.add_nodes_from(bound_vertices)
    graph.add_edges_from(bound_edges, real=1)
    # Return the added vertex labels and edges
    return bound_vertices, bound_edges


def get_face_sight_cells(selected_faces, ordered_face_edges, graph, positions,
                         bounds=((-1, -1), (-1, 1), (1, 1), (1, -1)),
                         outer=False):
    #
    all_face_edges = unlist([ordered_face_edges.get(face) for face in selected_faces])
    face_edge_map = {edge: [edge] for edge in all_face_edges}
    sight_cells = {face: None for face in selected_faces}

    # Iterate over all faces
    for face in selected_faces:

        # Get Vertices and ensure they are listed in counter-clockwise order
        face_edges = ordered_face_edges[face]
        face_vertices = get_sorted_face_vertices(face_edges, is_sorted=True)
        if calculate_face_signed_area(face_vertices, positions) < 0:
            face_vertices = list(reversed(face_vertices))

        # Calculate Inner Angles to check convexity
        face_angles = calculate_face_inner_angles(face_vertices, positions) if not outer \
            else calculate_face_outer_angles(face_vertices, positions)

        # If the face is convex, the face is the sight-cell
        if is_convex(face_angles) and not outer:
            face_edge_map[face].update({edge: [edge] for edge in face_edges})
            sight_cells[face] = face_edges

        # Otherwise, split into sight cells and return set of frozen sets of vertices per sight-cell
        else:

            # Project sight lines if possible
            edge_to_virtual_vertices, added_vertices = project_face_against_self(
                face, ordered_face_edges, face_edge_map, graph, positions, bounds, outer=False)

            # Update Graph and Edge Map
            cells = update_graph_and_virtual_edge_map(face, added_vertices, ordered_face_edges, face_edge_map,
                                                      edge_to_virtual_vertices, graph, positions, outer=False)

            # Log Identified Sight Cells
            sight_cells[face] = cells

    # Return dictionary of sight cells per face and a map of edges to virtual edges
    return sight_cells, face_edge_map


def get_clockwise_face_vertices(face, ordered_face_edges, face_edge_map, positions, original=False):
    if not original:
        face_edges = unlist([face_edge_map.get(edge) for edge in ordered_face_edges[face]])  # TODO: edges busted
    else:
        face_edges = [edge for edge in ordered_face_edges[face]]
    face_vertices = get_sorted_face_vertices(face_edges, is_sorted=True)
    if calculate_face_signed_area(face_vertices, positions) < 0:
        face_vertices = list(reversed(face_vertices))
    return face_vertices


def update_graph_and_virtual_edge_map(face, added_vertices, ordered_face_edges, face_edge_map,
                                      edge_to_virtual_vertices, graph, positions, outer=False):

    # Update List of Vertices with added vertices and the set of edges which define the face
    face_vertices = get_clockwise_face_vertices(face, ordered_face_edges, face_edge_map, positions)
    candidate_edges = unlist([face_edge_map.get(edge) for edge in face_edge_map.keys()])
    face_vertices = face_vertices + added_vertices + list(set(unlist(candidate_edges)))

    # Update the Graph and Its Positions
    cells, virtual_edge_map = update_sight_line_graph(edges=candidate_edges,
                                                      face_vertices=face_vertices,
                                                      edge_to_virtual_vertices=edge_to_virtual_vertices,
                                                      graph=graph,
                                                      positions=positions,
                                                      outer=outer)

    # Update the map of real to virtual edge maps
    deep_update_of_virtual_edge_map(complete_map=face_edge_map,
                                    new_map=virtual_edge_map)

    return cells


def project_face_against_self(face, ordered_face_edges, face_edge_map, graph, positions, bounds, outer=False):
    # Get the set of vertices which define the current face
    face_vertices = get_clockwise_face_vertices(face, ordered_face_edges, face_edge_map, positions, original=True)

    # Calculate Inner Angles to identify joint vertices
    face_angles = calculate_face_outer_angles(face_vertices, positions)

    # Replace original edges with their virtual counterparts
    candidate_edges = unlist([face_edge_map.get(edge) for edge in face_edge_map.keys()])

    # Extend all sight lines, add new vertices where necessary, and keep track of edge bisection
    edge_to_virtual_vertices, added_vertices = project_face_sight_lines(edges=candidate_edges,
                                                                        vertices=face_vertices,
                                                                        inner_angles=face_angles,
                                                                        edge_map=face_edge_map,
                                                                        graph=graph,
                                                                        positions=positions,
                                                                        bounds=bounds,
                                                                        outer=outer)

    # Return
    return edge_to_virtual_vertices, added_vertices


def project_outer_face_against_singleton():
    # extend_sight_line(joint_vertex, connecting_vertex, inner_angles, vertices, edges, graph, positions, bounds, outer)
    pass


def project_outer_face_against_non_face():
    pass


def project_outer_face_against_another_face(face, other_face, edge_map, ordered_face_edges, positions, graph, bounds):

    # Replace original edges with their virtual counterparts
    candidate_edges = unlist([edge_map.get(edge) for edge in edge_map.keys()])
    other_face_vertices = get_clockwise_face_vertices(
        other_face, ordered_face_edges, edge_map, positions, original=True)
    other_face_angles = calculate_face_outer_angles(other_face_vertices, positions)

    # Get Vertices and ensure they are listed in counter-clockwise order
    face_vertices = get_clockwise_face_vertices(face, ordered_face_edges, edge_map, positions, original=True)
    face_angles = calculate_face_outer_angles(face_vertices, positions)

    #
    edge_to_virtual_vertices, added_vertices = project_additional_sight_lines(edges=candidate_edges,
                                                                              origin_vertices=face_vertices,
                                                                              origin_angles=face_angles,
                                                                              target_vertices=other_face_vertices,
                                                                              edge_map=edge_map,
                                                                              target_angles=other_face_angles,
                                                                              graph=graph,
                                                                              positions=positions,
                                                                              bounds=bounds,
                                                                              outer=True)

    return edge_to_virtual_vertices, added_vertices


def get_subgraph(nodes, edges, graph, positions):
    sub_graph, sub_positions = copy.deepcopy(graph), copy.deepcopy(positions)
    sub_graph.remove_nodes_from([n for n in graph if n not in set(nodes)])
    [sub_positions.pop(n, None) for n in graph if n not in nodes]
    sub_graph.remove_edges_from([e for e in sub_graph.edges if set(e) not in [set(edge) for edge in edges]])
    return sub_graph, sub_positions


def find_outer_face_sight_cells(selected_faces, ordered_face_edges, graph, positions, is_cycle, bounds):
    # TODO: a bisected REAL edge will not be extended since we are looking up the original edge sets, whcih don't
    # TODO: exist anymore, i think. look at (7, 14) and (14, 8) not being extended in 1.5

    # Create lists of vertices and edges that define the outer face
    all_face_edges = unlist([ordered_face_edges.get(face) for face in selected_faces if len(face) > 1])
    outer_face_vertices = list(set(unlist(all_face_edges)))

    # Create A graph consisting only of outer-face defining vertices and edges
    outer_graph, outer_positions = get_subgraph(outer_face_vertices, all_face_edges, graph, positions)

    # Define the embedding boundary and the face edge map
    bound_vertices, bound_edges = add_boundary_to_graph(bounds, outer_graph, outer_positions)
    face_edge_map = {edge: [edge] for edge in all_face_edges + bound_edges}

    # Iterate over all faces
    selected_face_list = list(selected_faces)
    for face_index, face in enumerate(selected_face_list):
        if not is_cycle[face]:
            continue

        # Add additional candidate edges if we are dealing with the outer face
        other_faces = copy.copy(selected_faces)
        other_faces.remove(face)

        # Project sight-lines within the currently selected face against itself
        edge_to_virtual_vertices, added_vertices = project_face_against_self(
            face, ordered_face_edges, face_edge_map, outer_graph, outer_positions, bounds, outer=True)

        # Update Graph and Virtual Edge Map with New added vertices
        update_graph_and_virtual_edge_map(face, added_vertices, ordered_face_edges, face_edge_map,
                                          edge_to_virtual_vertices, outer_graph, outer_positions, outer=True)

        # DEBUG: Draw Initial Embedding
        draw_graph(graph=outer_graph, positions=outer_positions)
        save_drawn_graph(f"./graph_{face}.png")

        # Iterate over all other faces that form the outer face
        for other_face in selected_faces:
            if other_face == face:
                continue

            # Check if the Other face is a cycle or not -> check angle visibilities + edge crossing visibilities
            if is_cycle[other_face]:
                edge_to_virtual_vertices, added_vertices = project_outer_face_against_another_face(
                    face, other_face, face_edge_map, ordered_face_edges, outer_positions, outer_graph, bounds)
            elif not is_cycle[other_face] and len(other_face) > 1:
                # TODO: if other face is not a cycle -> special visiblility checking of ONLY edge intersections
                continue
            elif not is_cycle[other_face] and len(other_face) == 1:
                # TODO: if length(outer_face) == 1 -> only project single vertex
                continue
            else:
                sys.exit("Non Accounted for Constellation of face!")

            # Update Graph and Virtual Edge Map with New added vertices
            update_graph_and_virtual_edge_map(face, added_vertices, ordered_face_edges, face_edge_map,
                                              edge_to_virtual_vertices, outer_graph, outer_positions, outer=True)

            # Draw Initial Embedding
            draw_graph(graph=outer_graph, positions=outer_positions)
            save_drawn_graph(f"./graph_{face}+{other_face}.png")

    # Identify all faces (i.e. sight cells in outer face)
    face_graph = nx.Graph(outer_graph.subgraph(nodes=list(set(unlist(unlist(list(face_edge_map.values())))))))
    sight_cells = find_inner_faces(face_graph, positions=outer_positions)

    # Remove the original set of faces that defined the outer face
    [sight_cells.remove(cell) for cell in copy.copy(sight_cells) for face in selected_faces if face.issubset(cell)]

    # Return the identified sight cells and the subgraph
    return sight_cells, face_edge_map, outer_graph, outer_positions


def deep_update_of_virtual_edge_map(complete_map, new_map):

    # TODO: edges must be frozensets -> mapping of (x, y) does not work for (y, x)
    #  maybe revamp this whole framework to just work on frozensets instead of tuples?
    for virtual_edge in new_map.keys():
        is_mapped = lambda edge, edge_list: frozenset(list(edge)) in [frozenset(e) for e in edge_list]
        mapped = [real_edge for real_edge in complete_map.keys() if is_mapped(virtual_edge, complete_map[real_edge])]

        if len(mapped) == 1:
            real_edge = mapped[0]
            real_edge_set_map = [frozenset(edge) for edge in complete_map[real_edge]]
            index = real_edge_set_map.index(frozenset(virtual_edge))
            real_edge_set_map[index:index] = new_map[virtual_edge]
            real_edge_set_map.remove(frozenset(virtual_edge))
            complete_map[real_edge] = [tuple(edge) for edge in real_edge_set_map]
        elif len(mapped) > 1:
            sys.exit("SHIT'S FUCKED")
        else:
            complete_map.update({virtual_edge: new_map[virtual_edge]})


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


def are_vertices_adjacent_virtually(vertex_a, vertex_b, edge_map):
    return ((vertex_a, vertex_b) in edge_map.keys()) or ((vertex_b, vertex_a) in edge_map.keys())
