import numpy as np


def calculate_face_centroid(face_vertex_positions):
    return sum(face_vertex_positions) / len(face_vertex_positions)


def get_split_vertex_locations(positions, target_face_subfaces):
    centroids = {target_face: dict() for target_face in target_face_subfaces.keys()}
    for target_face in target_face_subfaces.keys():
        for subface_set in target_face_subfaces[target_face]:
            subface_positions = [positions[vertex] for vertex in list(subface_set)]
            centroids[target_face][subface_set] = calculate_face_centroid(subface_positions)
    return centroids

