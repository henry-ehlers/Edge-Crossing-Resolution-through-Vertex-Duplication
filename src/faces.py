from src.edges import *

import matplotlib.path as mpltPath
import networkx as nx
import itertools as it
import numpy as np
import copy
import sys


def find_all_subfaces(graph, virtual_edge_set_map, target_face_to_vertex_map):

    # Prepare dictionary of sets within which to store all found faces per face
    subfaces = {target_face: set() for target_face in target_face_to_vertex_map.keys()}

    # Limit search of cycle basis to subfaces
    for target_face in target_face_to_vertex_map.keys():
        vertices_to_keep = set()

        # Mark vertices based on their inclusion in virtual edge sets pertaining to the current set
        for edge_set in virtual_edge_set_map:
            intersection = target_face_to_vertex_map[target_face] & edge_set
            if len(intersection) >= 2:
                [vertices_to_keep.add(vertex) for vertex in list(edge_set)]
            elif len(intersection) == 1:
                vertices_to_keep.add(intersection.pop())
            else:
                continue

        # Keep only vertices pertaining to current subgraph and search for cycle basis within
        subgraph = copy.deepcopy(graph).subgraph(list(vertices_to_keep))
        subfaces[target_face] = find_all_faces(subgraph)  # todo: also pass positions?

    # Return all subfaces for each target face
    return subfaces


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


def get_ordered_face_edges(faces, graph):
    ordered_face_edges = dict.fromkeys(faces)
    for face in faces:
        ordered_face_edges[face] = get_face_vertex_sequence(face, graph)
    return ordered_face_edges


def shrink_cycle(cycle, other_cycles, sorted_edges, graph, positions):

    sorted_vertices = get_sorted_face_vertices(sorted_edges[cycle], is_sorted=True)
    cycle_coordinates = [positions[vertex] for vertex in sorted_vertices]

    # Get positions from polygon of ordered points of current cycle
    cycle_path = mpltPath.Path(cycle_coordinates[0:-1])
    for other_cycle in other_cycles:

        # Skip if the cycle == the other cycle
        if other_cycle == cycle:
            continue

        vertex_intersection = cycle.intersection(other_cycle)
        if len(vertex_intersection) == 0:
            continue

        remaining_vertices = other_cycle - vertex_intersection
        if len(remaining_vertices) == 0:
            continue

        remaining_coordinates = [positions[vertex] for vertex in remaining_vertices]
        in_side = cycle_path.contains_points(remaining_coordinates)
        if not all(in_side):
            if not all(not in_side[index] for index in range(0, len(remaining_vertices))):
                print("SHRINKING CYCLES are fucked")
            continue

        #
        new_cycle, new_edge_list = None, None
        if len(vertex_intersection) == 1:
            # TODO: fix this section -> edge order is unclear
            new_cycle = cycle.union(remaining_vertices)
            new_edge_list = get_face_vertex_sequence(new_cycle, graph)

        elif len(vertex_intersection) >= 2:
            cycle_edge_list = set([frozenset(edge) for edge in sorted_edges[cycle]])
            other_edge_list = set([frozenset(edge) for edge in sorted_edges[other_cycle]])

            new_edge_set = cycle_edge_list.symmetric_difference(other_edge_list)
            new_edge_list = sort_face_edges(list([tuple(edge) for edge in new_edge_set]))
            new_cycle = frozenset().union(*new_edge_set)

        # Mention that something is broken
        if (new_cycle is None) or (new_edge_list is None):
            sys.exit("CYCLES FUCKED")

        # Update the list of cycles
        other_cycles.remove(cycle)
        other_cycles.add(new_cycle)

        # Update the cycle sorted edge dictionary
        sorted_edges[new_cycle] = new_edge_list
        sorted_edges.pop(cycle)

        # Recurse and return
        return shrink_cycle(new_cycle, other_cycles, sorted_edges, graph, positions)


def find_singleton_cycles(cycles, graph, as_set=True):

    # Get Ordered Face edges for all identified cycles
    ordered_face_edges = get_ordered_face_edges(cycles, graph)

    # Find edges not Featured in Cycles
    non_cycle_edges = set()
    for edge in list(graph.edges):
        in_cycle = [set(edge) in [set(cycle_edge) for cycle_edge in ordered_face_edges[cycle]] for cycle in cycles]
        if any(in_cycle):
            continue
        non_cycle_edges.add(frozenset(edge))

    # Add found singleton cycles to identified ones
    if as_set:
        [cycles.add(singleton_cycle) for singleton_cycle in non_cycle_edges]
    else:
        [cycles.append(singleton_cycle) for singleton_cycle in non_cycle_edges]


def find_all_faces(graph, positions=None, as_set=True):

    # TODO: in all subsequent functions which use the face dictionary, check the length of the face (no singletons)

    # Identify the minimum cycle basis of the graph
    cycles = nx.minimum_cycle_basis(G=graph)
    cycles = set([frozenset(cycle) for cycle in cycles]) if as_set else [set(cycle) for cycle in cycles]

    # Check for and add Singleton edges as singleton cycles
    find_singleton_cycles(cycles, graph, as_set)

    # If no positions are given with which to check the validity of the cycles, return the cycles as faces
    if positions is None:
        return cycles

    # Get the ordered edge list per cycle
    cycle_ordered_edges = get_ordered_face_edges(cycles, graph)
    # For each Cycle recursively check that the face is indeed a face
    faces = copy.copy(cycles)
    for cycle in cycles:
        shrink_cycle(cycle=cycle,
                     other_cycles=faces,
                     sorted_edges=cycle_ordered_edges,
                     graph=graph,
                     positions=positions)

    # Return set of faces (frozen sets of vertices)
    return faces


def color_selected_faces(graph, face_set, face_edge_map):
    for face in face_set:
        for edge in face_edge_map[frozenset(face)]:
            graph.edges[edge]["target"] = 1
    return graph


def close_closed_face(convex_face_edge_list, graph):

    # Note down the face's edges
    face_edges = []

    # Iterate over all ORDERED vertices that make up the new convex face
    for index in range(0, len(convex_face_edge_list)):

        # Extract two adjacent vertices from the face
        vertex_a = convex_face_edge_list[index]
        vertex_b = convex_face_edge_list[(index + 1) % len(convex_face_edge_list)]

        # Check if an edge already exists between these vertices
        try:
            edge = graph.edges[vertex_a, vertex_b]
        except KeyError:
            edge = None

        # If no edge exists, add one and update the convex face's ordered edge list
        if edge is None:
            graph.add_edge(vertex_a, vertex_b, virtual=1, segment=0, target=0)

        # Append (new) edge to this face's list of edges
        face_edges.append((vertex_a, vertex_b))

    # Return update list of unsorted face edges
    return face_edges


def find_outer_face(ordered_face_edges, graph):

    # Initialize a dictionary which maps edges to the number of faces they are in
    faces_per_edge = dict.fromkeys(list(graph.edges()), 0)

    # Iterate over all faces and increment counts of edges within them
    for face in ordered_face_edges.keys():
        for edge in ordered_face_edges[face]:
            edge_a, edge_b = (edge[0], edge[1]), (edge[1], edge[0])
            if edge_a in faces_per_edge.keys():
                faces_per_edge[edge_a] += 1
            elif edge_b in faces_per_edge.keys():
                faces_per_edge[edge_b] += 1
            else:
                sys.exit("Shit's fucked")

    # Find edges which only map to a single face
    # min_face_count = min(faces_per_edge.values())
    faces = set([frozenset(edge) for edge in faces_per_edge.keys() if faces_per_edge[edge] == 1])

    find_vertex_sets_from_edges(faces)

    # Return unique vertex sets from the found singleton edges
    return faces


def find_vertex_sets_from_edges(edge_sets):
    for edge_set_a in edge_sets:
        for edge_set_b in edge_sets:
            if edge_set_a == edge_set_b:
                continue
            if not edge_set_a.intersection(edge_set_b):
                continue

            new_edge_set = edge_set_a.union(edge_set_b)
            edge_sets.remove(edge_set_a)
            edge_sets.remove(edge_set_b)
            edge_sets.add(new_edge_set)

            return find_vertex_sets_from_edges(edge_sets)


def build_face_to_edge_map(graph, faces):

    # TODO: the face set now consists of shortest cycles, so this basic look-up should be acceptable for now
    # TODO: however, for outer faces, there may be cycles within the the defined cycle that should not be included
    face_edge_map = dict()
    for face in faces:
        face_edge_map[face] = list()
        face_list = list(it.combinations(list(face), 2))
        for vertex_combination in face_list:
            vertex_a = vertex_combination[0]
            vertex_b = vertex_combination[1]
            found_edge = get_edge_connecting_vertices(graph, vertex_a, vertex_b)
            if found_edge: face_edge_map[face].append(found_edge)

    return face_edge_map


def get_edge_connecting_vertices(graph, vertex_a, vertex_b):
    for edge in graph.edges:
        if (vertex_a == edge[0] and vertex_b == edge[1]) or (vertex_b == edge[0] and vertex_a == edge[1]):
            return edge


def find_face_vertex_incidence(faces, target_vertices):

    # Initialize an empty dictionary (using sets as keys) to store vertex incidence sets
    face_incidences = dict()  # face_incidences = {face: dict() for face in faces}

    # Initialize set of vertex set as list, and target vertices as set
    faces = list([list(face) for face in faces])
    target_vertex_set = set(target_vertices)

    # Iterate over all faces
    for face_index_a in range(0, len(faces) - 1):

        # Extract current face, transform back into set
        face_a = frozenset(faces[face_index_a])
        face_incidences[face_a] = dict()

        # Determine how many target vertices are incident and left over
        incidence_a = face_a & target_vertex_set
        remaining_targets = target_vertex_set - incidence_a

        # Iterate over remaining faces
        for face_index_b in range(face_index_a + 1, len(faces)):

            # Extract current face, transform back into set
            face_b = frozenset(faces[face_index_b])

            # Determine how many remaining target faces are incident
            incidence_b = face_b & remaining_targets
            incident_vertices = incidence_a.union(incidence_b)

            # Store results only if non-empty
            if incident_vertices:
                face_incidences[face_a][face_b] = incident_vertices

    # Return Dictionary of face vertex incidence
    return face_incidences


def get_maximally_incident_faces(face_incidences):
    max_incidence = float("-inf")
    selected_faces = []
    for face_a in face_incidences.keys():
        for face_b in face_incidences[face_a].keys():
            vertex_incidence = len(face_incidences[face_a][face_b])
            if vertex_incidence > max_incidence:
                max_incidence = vertex_incidence
                selected_faces = [[face_a, face_b]]
            elif vertex_incidence == max_incidence:
                selected_faces.append([face_a, face_b])

    return max_incidence, selected_faces


def cross_product(vector_a, vector_b):
    return (vector_a[0] * vector_b[1]) - (vector_b[0] * vector_a[1])


def vector_angle(vector_1, vector_2):

    # Calculate Unit Vectors of Input Vectors
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)

    # Calculate Dot Product and Signed Angle in Radians
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(dot_product)

    # Return Angle in Degrees
    return np.degrees(angle)


# Calculate Signed Area of Polygon
def signed_area(ordered_points_list):
    x, y = [point[0] for point in ordered_points_list], [point[1] for point in ordered_points_list]
    return sum(x[i] * (y[i + 1] - y[i - 1]) for i in range(-1, len(x) - 1)) / 2.0


def calculate_face_signed_area(ordered_face_vertices, positions):
    vertex_points = [positions[vertex] for vertex in ordered_face_vertices]
    area = signed_area(vertex_points)
    return area


def calculate_outer_angle(point_a, point_b, point_c):

    # Assumed that point b connects to both a and c
    vector_1, vector_2 = point_a - point_b, point_c - point_b

    # Calculate Signed Angle between two Vectors
    signed_angle = vector_angle(vector_1, vector_2)

    # Return Outer angle
    return signed_angle if cross_product(vector_1, vector_2) < 0.0 else 360 - signed_angle


def calculate_face_outer_angles(counter_clockwise_face_vertices, positions):

    # Initialize emtpy dictionary to store inner angles
    outer_angles = {vertex: None for vertex in counter_clockwise_face_vertices}

    # Iterate over all vertices, and calculate angle for each
    for vertex_index in range(0, len(counter_clockwise_face_vertices)):

        # Get three vertices that form an angle (the ordered hereof is crucial)
        vertex_a = counter_clockwise_face_vertices[vertex_index]
        vertex_b = counter_clockwise_face_vertices[vertex_index - 1]
        vertex_c = counter_clockwise_face_vertices[vertex_index - 2]

        # Calculate Inner angle and store as with center vertex as key
        point_a, point_b, point_c = positions[vertex_a], positions[vertex_b], positions[vertex_c]
        outer_angles[vertex_b] = calculate_outer_angle(point_a, point_b, point_c)

    # Return all inner angles
    return outer_angles


def calculate_inner_angle(point_a, point_b, point_c):

    # Assumed that point b connects to both a and c
    vector_1, vector_2 = point_a - point_b, point_c - point_b

    # Calculate Signed Angle between two Vectors
    signed_angle = vector_angle(vector_1, vector_2)

    # Return inner angle
    return signed_angle if cross_product(vector_1, vector_2) > 0.0 else 360 - signed_angle


def calculate_face_inner_angles(counter_clockwise_face_vertices, positions):

    # Initialize emtpy dictionary to store inner angles
    inner_angles = {vertex: None for vertex in counter_clockwise_face_vertices}

    # Iterate over all vertices, and calculate angle for each
    for vertex_index in range(0, len(counter_clockwise_face_vertices)):

        # Get three vertices that form an angle (the ordered hereof is crucial)
        vertex_a = counter_clockwise_face_vertices[vertex_index]
        vertex_b = counter_clockwise_face_vertices[vertex_index - 1]
        vertex_c = counter_clockwise_face_vertices[vertex_index - 2]

        # Calculate Inner angle and store as with center vertex as key
        point_a, point_b, point_c = positions[vertex_a], positions[vertex_b], positions[vertex_c]
        inner_angles[vertex_b] = calculate_inner_angle(point_a, point_b, point_c)

    # Return all inner angles
    return inner_angles

