from py2d.Math import Polygon, Vector
from src.edges import *
import networkx as nx
import itertools as it
import numpy as np
import copy


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


def is_face_convex(face, ordered_face_edges, position):
    if len(face) == 3:
        return True, None
    else:
        position_list = [tuple(position[edge[0]]) for edge in ordered_face_edges]
        face_polygon = Polygon.from_tuples(position_list)
        if face_polygon.is_convex():
            return True, None
        else:
            return False, Polygon.convex_decompose(face_polygon)


def make_faces_convex(faces, ordered_face_edges, graph, positions):
    convex_faces = set()
    for face in faces:
        edge_list = ordered_face_edges[face]
        is_convex, decompositions = is_face_convex(face, edge_list, positions)
        if is_convex:
            convex_faces.add(face)
        else:
            coordinate_decompositions, index_decompositions = decompositions
            print(f"non_convex_face edge list: {edge_list}")
            decomposition = [[]] * len(index_decompositions)
            for index, index_decomposition in enumerate(index_decompositions):
                vertex_decomposition = [edge_list[index][0] for index in index_decomposition]
                convex_faces.add(frozenset(vertex_decomposition))
                decomposition[index] = vertex_decomposition
                close_closed_face(vertex_decomposition, graph)
            [print(decomp) for decomp in decomposition]

    # TODO: add edges where needed between decompositions
    return convex_faces


def close_closed_face(convex_face_edge_list, graph):

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

        # If no edge exists, add one
        if edge is None:
            graph.add_edge(vertex_a, vertex_b, virtual=1, segment=0, target=0)


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
