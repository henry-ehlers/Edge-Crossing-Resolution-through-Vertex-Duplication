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

    faces = list([list(face) for face in faces])
    target_vertex_set = set(target_vertices)
    print(len(faces))

    for face_index_a in range(0, len(faces) - 1):
        face_a = faces[face_index_a]
        face_incidences[frozenset(face_a)] = dict()
        print("\n>{} - {}".format(face_index_a, range(face_index_a+1, len(faces))))
        incidence_a = set(face_a) & target_vertex_set
        remaining_targets = target_vertex_set - incidence_a
        # print("-------------------------------------------")
        # print("face a:      {}".format(face_a))
        # print("target:      {}".format(target_vertex_set))
        # print("incidence_a: {}".format(incidence_a))
        # print("remaining:   {}".format(remaining_targets))
        for face_index_b in range(face_index_a + 1, len(faces)):
            face_b = faces[face_index_b]
            incidence_b = set(face_b) & remaining_targets
            # print("face b:      {}".format(face_b))
            # print("target       {}".format(remaining_targets))
            # print("incidence_b: {}".format(incidence_b))
            # print("remaining:   {}".format(remaining_targets - incidence_b))
            incident_vertices = incidence_a.union(incidence_b)
            if incident_vertices:
                face_incidences[frozenset(face_a)][frozenset(face_a)] = incident_vertices
                print("final:       {}".format(face_incidences[frozenset(face_a)][frozenset(face_a)]))
        print(len(faces))
    return face_incidences
