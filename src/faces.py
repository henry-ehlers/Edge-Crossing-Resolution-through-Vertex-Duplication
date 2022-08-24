from src.edges import *
from src.line_segments import *
from src.edge_crossings import *
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
        subfaces[target_face] = find_all_faces(subgraph)

    # Return all subfaces for each target face
    return subfaces


def get_ordered_face_edges(faces, graph):
    ordered_face_edges = dict.fromkeys(faces)
    for face in faces:
        ordered_face_edges[face] = get_face_vertex_sequence(face, graph)
    return ordered_face_edges


def find_all_faces(graph):

    # Identify the minimum cycle basis of the graph
    # TODO: this does not identify outer faces (i.e. the outer cycle)
    cycles = nx.minimum_cycle_basis(G=graph)

    # Convert the minimum cycle basis into face sets
    # TODO: check that cycles are indeed legal faces
    faces = set([frozenset(cycle) for cycle in cycles])

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
            print(f"Added {vertex_a} and {vertex_b}")
            graph.add_edge(vertex_a, vertex_b, virtual=1, segment=0, target=0)

        # Append (new) edge to this face's list of edges
        face_edges.append((vertex_a, vertex_b))

    # Return update list of unsorted face edges
    return face_edges


def find_outer_face(ordered_face_edges, graph):
    faces_per_edge = dict.fromkeys(list(graph.edges()), 0)
    for face in ordered_face_edges.keys():
        for edge in ordered_face_edges[face]:
            edge_a, edge_b = (edge[0], edge[1]), (edge[1], edge[0])
            if edge_a in faces_per_edge.keys():
                faces_per_edge[edge_a] += 1
            elif edge_b in faces_per_edge.keys():
                faces_per_edge[edge_b] += 1
            else:
                sys.exit("Shit's fucked")

    min_face_count = min(faces_per_edge.values())
    print(f"minimum face count: {min_face_count}")
    outer_face = [edge for edge in faces_per_edge.keys() if faces_per_edge[edge] == min_face_count]
    return outer_face


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


def get_sight_cells(faces, ordered_face_edges, graph, positions, bounds=((-1, -1), (-1, 1), (1, 1), (1, -1))):

    #
    sight_cells = {face: None for face in faces}

    # Iterate over all faces
    for face in faces:

        # Get Vertices and ensure they are listed in counter-clockwise order
        face_edges = ordered_face_edges[face]
        face_vertices = get_sorted_face_vertices(face_edges, is_sorted=True)
        if calculate_face_signed_area(face_vertices, positions) < 0:
            face_vertices = list(reversed(face_vertices))

        # Calculate Inner Angles
        face_inner_angles = calculate_face_inner_angles(face_vertices, positions)
        if is_convex(face_inner_angles):
            sight_cells[face] = face_edges
        else:
            # split into sight cells
            # TODO: can we recycle the calculated inner_angles somehow? we need to know the reference and both endpoints
            split_face_into_sight_cells(face_edges, face_vertices, graph, positions)

    return sight_cells


def is_convex(inner_angles):
    return all(angle <= 180.0 for angle in inner_angles.values())


def split_face_into_sight_cells(edges, vertices, graph, positions, bounds=((-1, -1), (-1, 1), (1, 1), (1, -1))):

    print(vertices)
    for origin_index in range(0, len(vertices)):

        for target_index in range(0, len(vertices)):

            vertex_a, vertex_b = vertices[origin_index], vertices[target_index]
            vertex_set = {vertex_a, vertex_b}

            print(f"\nVertex Set: {vertex_set}")
            if (vertex_a == vertex_b): continue

            ind_diff = abs(origin_index - target_index)
            adjacent = True if ((ind_diff == 1) or (ind_diff == len(vertices) - 1)) else False
            is_visible = is_vertex_visible(vertex_a, vertex_b, edges, positions)
            print(f"Vertices {vertex_a} and {vertex_b} can see one-another: {is_visible}")

            # TODO: if is visible, extend line beyond until it intersect another line
            # TODO: Check using angle visibility in which direction to extend
            # TODO: if multiple intersections, take closest one (the other is not visible)

            # TODO: MUST CONSIDER EXTENSIONS IN BOTH DIRECTIONS TO COVER EDGE_CASES, REGARDLESS OF WHETHER ADJACENT OR NOT
            if is_visible:
                extend_sight_line(vertex_a, vertex_b, edges, graph, positions, bounds)


def squared_distance(point_a, point_b):
    (x1, y1) = point_a
    (x2, y2) = point_b
    return (x2 - x1) ** 2 + (y2 - y1) ** 2


def extend_sight_line(vertex_a, vertex_b, adjacent, edges, graph, positions, bounds):

    # Calculate intersections of extended line with boundaries in both directions
    bound_intersections = extend_line(positions[vertex_a], positions[vertex_b], bounds)
    print(f"intersections: {bound_intersections}")

    # If the vertices are not adjacent, check which one is legal based on
    if not adjacent:
        intersection_reference = bound_intersections[0] \
            if check_vertex_visibility_by_angle(vertex_a, bound_intersections[0], edges, positions) \
            else bound_intersections[1]
    else:
        distance_from_a = [squared_distance(vertex_a, intersection) for intersection in bound_intersections]
        distance_from_b = [squared_distance(vertex_b, intersection) for intersection in bound_intersections]
        intersection_reference = bound_intersections[0] \
            if (distance_from_b[0] < distance_from_b[1]) and (distance_from_a[0] > distance_from_a[1]) \
            else bound_intersections[1]
    print(f"reference: {intersection_reference}")
    possible_crossing_edges = [edge for edge in edges if (vertex_a not in edge) and (vertex_b not in edge)]
    print(f"possible crossing edges: {possible_crossing_edges}")
    point_a, point_b = positions[vertex_a], intersection_reference
    print(f"Points A and B: {point_a} and {point_b}")
    for edge in possible_crossing_edges:
        point_c, point_d = positions[edge[0]], positions[edge[1]]
        print(f"Points C and D: {point_c} and {point_d}")
        edge_intersections = line_intersection(point_a, point_b, point_c, point_d)
        if edge_intersections is not None:
            print("INTERSECTION")
    pass


def is_vertex_visible(vertex_a, vertex_b, edges, positions):

    # Check Angle first, and only if the angle is legal, exclude possible intersection edges
    angle_visibility = check_vertex_visibility_by_angle(vertex_a, vertex_b, edges, positions)
    if not angle_visibility:
        return angle_visibility

    # Check specific edge crossings between (Vertex A, Vertex B) and all non-incident edges
    possible_crossing_edges = [edge for edge in edges if (vertex_a not in edge) and (vertex_b not in edge)]
    crossing_visibility = check_vertex_visibility_by_crossing(vertex_a, vertex_b, possible_crossing_edges, positions)
    if not crossing_visibility:
        return crossing_visibility

    # All Conditions met for Visibility
    return True


def check_vertex_visibility_by_angle(vertex_a, vertex_or_position_b, edges, positions):

    # Calculate Angle of between Vertex A and its two neighbors relative to Vertex Edge_a[0]
    edge_a = [edge for edge in edges if vertex_a == edge[1]][0]  # Order of vertices in edges important; here 2nd
    edge_b = [edge for edge in edges if vertex_a == edge[0]][0]  # Order of vertices in edges important; here 1st
    point_a, point_b, point_c = positions[edge_a[0]], positions[edge_a[1]], positions[edge_b[1]]
    observed_angle = calculate_inner_angle(point_a, point_b, point_c)
    # TODO: make sure the observed angle actually IS the inner angle; before we needed to use the signed area

    # Calculate angle between Vertex A and Vertex B relative to Vertex Edge_a[0]
    point_a, point_b = positions[edge_a[0]], positions[edge_a[1]],
    point_c = positions[vertex_or_position_b] if isinstance(vertex_or_position_b, int) else vertex_or_position_b
    hypothetical_angle = calculate_inner_angle(point_a, point_b, point_c)

    # If the angle between Vertex A and B is larger than between Vertex A and its neighbors, Vertex B is not visible
    return False if hypothetical_angle > observed_angle else True


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
