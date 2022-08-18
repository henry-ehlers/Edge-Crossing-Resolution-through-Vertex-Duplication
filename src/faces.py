from py2d.Math import Polygon, Vector
from src.edges import *
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

    # Initialize new containers for convex faces and their sorted edges
    convex_faces, convex_ordered_face_edges = set(), dict()

    # Iterate over all faces to check whether they are convex
    for face in faces:

        # Extract ordered edge list and (with vertex positions) check convexity
        edge_list = ordered_face_edges[face]
        is_convex, decompositions = is_face_convex(face, edge_list, positions)

        # If the face is convex, copy content
        if is_convex:
            convex_ordered_face_edges[face] = ordered_face_edges[face]
            convex_faces.add(face)

        # If face is not convex, decompose it into distinct convex polygonal faces
        else:
            coordinate_decompositions, index_decompositions = decompositions
            print(f"index_decompositions: {index_decompositions}")

            # Iterate over all identified convex faces
            for index, index_decomposition in enumerate(index_decompositions):

                # Get edges and vertices of current convex face
                old_sorted_edges = [edge_list[edge_index] for edge_index in index_decomposition]
                vertex_decomposition = [face_edge[0] for face_edge in old_sorted_edges]

                # Add virtual edges to close the created convex faces
                convex_face_edges = close_closed_face(vertex_decomposition, graph)
                print(f"edge: {convex_face_edges}")

                # Sort edges and add to convex face sorted edge dictionary
                sorted_convex_face_edges = sort_face_edges(convex_face_edges)
                print(f"sorted convex vertices: {sorted_convex_face_edges}")

                # Add vertices and edges to convex sets
                convex_ordered_face_edges[frozenset(vertex_decomposition)] = sorted_convex_face_edges
                convex_faces.add(frozenset(vertex_decomposition))

    # Return new containers of convex faces and their sorted edges
    return convex_faces, convex_ordered_face_edges


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
