from src.edge_crossings import *
import pandas as pd


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


def select_best_vertex_split_pair(induced_edge_crossings):

    
    for target_face_a in induced_edge_crossings.keys():
        for target_face_b in induced_edge_crossings.keys():
            if target_face_a == target_face_b:
                continue
    return None