from src.edge_crossings import *
import pandas as pd
import itertools


def calculate_face_centroid(face_vertex_positions):
    return sum(face_vertex_positions) / len(face_vertex_positions)


def get_split_vertex_locations(positions, target_face_subfaces):
    centroids = {target_face: dict() for target_face in target_face_subfaces.keys()}
    for target_face in target_face_subfaces.keys():
        for subface_set in target_face_subfaces[target_face]:
            subface_positions = [positions[vertex] for vertex in list(subface_set)]
            centroids[target_face][subface_set] = calculate_face_centroid(subface_positions)
    return centroids


def calculate_induced_edge_crossings(graph, positions, centroids, target_neighbors):

    # Extract all real edges
    real_edges = [edge for edge in graph.edges if graph.edges[edge]["virtual"] == 0]
    print(real_edges)

    # Store all edge crossings
    induced_edge_crossings = {target_face: dict() for target_face in centroids.keys()}

    # Iterate over all target faces
    for target_face in centroids.keys():

        # Iterate over all subfaces in the current target face
        for subface in centroids[target_face].keys():
            point_a = centroids[target_face][subface]

            # Store the number of intersections per neighbor
            intersections = [0] * len(target_neighbors)

            # Iterate over all neighbors that still need to be connected
            for neighbor_index, neighbor in enumerate(target_neighbors):
                point_b = positions[neighbor]

                # Iterate over all edges in the graph
                for edge in real_edges:

                    # Extract vertices and positions that form edge
                    vertex_c, vertex_d = edge[0], edge[1]
                    point_c, point_d = positions[vertex_c], positions[vertex_d]

                    # Calculate intersection and store if they do
                    if line_intersection(point_a, point_b, point_c, point_d) is not None:
                        intersections[neighbor_index] += 1

            # Store the number of edge crossings
            induced_edge_crossings[target_face][subface] = tuple(intersections)

    # Return Complete Set
    return induced_edge_crossings


def get_split_vertex_pairs(induced_edge_crossings):

    # Initialize dictionaries to store results
    pair_induced_crossings = dict()
    neighbor_assignment = dict()

    # Define all possible combinations of target faces and iterate over them
    target_face_pairs = itertools.combinations(induced_edge_crossings.keys(), 2)
    for face_pair in target_face_pairs:
        face_a, face_b = face_pair
        print(f"Faces {face_a} and {face_b}")

        pair_induced_crossings.update({subface: dict() for subface in induced_edge_crossings[face_a].keys()})
        neighbor_assignment.update({subface: dict() for subface in induced_edge_crossings[face_a].keys()})

        # Iterate over each target face's subfaces and extract their induced edge crossings numbers
        for subface_a in induced_edge_crossings[face_a].keys():
            crossings_a = induced_edge_crossings[face_a][subface_a]
            print(f"Subface-A: {subface_a}")

            for subface_b in induced_edge_crossings[face_b].keys():
                crossings_b = induced_edge_crossings[face_b][subface_b]
                print(f"Subface-B: {subface_b}")

                # Select the best arrangement of edge crossings
                best_crossings = [None] * len(crossings_a)
                assignment = [None] * len(crossings_a)
                for neighbor in range(0, len(crossings_a)):
                    # TODO: Resolve ties, maybe using edge length or distance from target face's centroid (FUNCTION)
                    neighbor_a_better = crossings_a[neighbor] <= crossings_b[neighbor]
                    best_crossings[neighbor] = crossings_a[neighbor] if neighbor_a_better else crossings_a[neighbor]
                    assignment[neighbor] = subface_a if neighbor_a_better else subface_b

                # Store results
                pair_induced_crossings[subface_a][subface_b] = sum(best_crossings)
                neighbor_assignment[subface_a][subface_b] = assignment

    return pair_induced_crossings, neighbor_assignment
