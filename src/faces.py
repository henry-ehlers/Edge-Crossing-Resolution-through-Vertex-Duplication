import networkx as nx
import numpy as np
import copy


def find_all_faces(graph):

    # Ensure the provided graph is plana
    is_planar, plane_graph = nx.check_planarity(graph)
    assert is_planar, \
        "Provided Graph is not planar."

    # Initialize set of all faces (sets of vertices)
    faces = set()

    # Iterate over each vertex in the drawing
    for origin_vertex in range(0, plane_graph.number_of_nodes()):

        # Collect the neighbors [w] of a vertex [v] in clock-wise order
        cw_neighbors = list(plane_graph.neighbors_cw_order(origin_vertex))

        # For each clockwise-sorted neighbor w, find the face to the right of the half-edge (v,w)
        for adjacent_vertex in cw_neighbors:
            face = frozenset(plane_graph.traverse_face(origin_vertex, adjacent_vertex))
            faces.add(face)

    # Return set of faces (frozen sets of vertices)
    return faces


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
