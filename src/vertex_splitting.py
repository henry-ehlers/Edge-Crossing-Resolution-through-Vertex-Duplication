from src.edge_crossings import *

import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import pandas as pd
import numpy as np
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
    real_edges = [edge for edge, real in nx.get_edge_attributes(G=graph, name="real").items() if real == 1]

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
                for (vertex_c, vertex_d) in real_edges:

                    # Do not check intersections with edges bounded by neighbor (ends count as crossings)
                    if neighbor in (vertex_c, vertex_d):
                        continue

                    # Extract vertices and positions that form edge
                    point_c, point_d = positions[vertex_c], positions[vertex_d]

                    # Calculate intersection and store if they do
                    if line_intersection(point_a, point_b, point_c, point_d) is not None:
                        intersections[neighbor_index] += 1

            # Store the number of edge crossings
            induced_edge_crossings[target_face][subface] = tuple(intersections)

    # Return Complete Set
    return induced_edge_crossings


def get_edge_crossing_table(induced_edge_crossings, target_neighbors):
    """
    a function to convert a provided dictionary of induced edge crossings to a dictionary of pandas dataframes.

    :param induced_edge_crossings: a nested dictionary which contains the number of edge crossings induced to connect
    the centroid of each face to each of the p target neighbors. The dictionary has a structure of
    {frozenset(face): frozenset(subface): (integer tuple of length p)}
    :param target_neighbors: a list of length p containing the vertices to which all centroids must connect

    :return: a dictionary of length f (where f is the number of faces), where each value is an n*p pandas dataframe,
    where n = the number of subfaces in face f, and p = the number of target neighbor vertices
    """

    # Initialize dictionary with two keys; one for each target face
    crossing_tables = {face: None for face in induced_edge_crossings.keys()}

    # Iter ate over the two target faces
    for target_face in crossing_tables.keys():

        # Get all induced edge crossings for all subfaces in current face
        crossing_tables[target_face] = get_subface_edge_crossing_table(
            face=target_face,
            subface_crossings=induced_edge_crossings[target_face],
            target_vertices=target_neighbors)

    # Return Edge Crossing Tables
    return crossing_tables


def get_subface_edge_crossing_table(face, subface_crossings, target_vertices):

    # Create Lists of sub-faces and faces
    sub_faces = [subface for subface in subface_crossings.keys()]
    faces = [face] * len(subface_crossings)

    # Tabulate the number of crossings for all target vertex and subface combinations
    crossings = np.empty(shape=(len(subface_crossings), len(target_vertices)),
                         dtype=int)
    for index, subface in enumerate(sub_faces):
        crossings[index, ...] = subface_crossings[subface]

    # Create crossing pandas data frame
    crossing_table = pd.DataFrame(data=crossings,columns=target_vertices, dtype=int)
    crossing_table.insert(loc=0, column="sub_face", value=sub_faces)
    crossing_table.insert(loc=0, column="face", value=faces)

    # Return table
    return crossing_table



def get_split_vertex_pairs(induced_edge_crossings):

    # Initialize dictionaries to store results
    pair_induced_crossings = dict()
    neighbor_assignment = dict()

    # Define all possible combinations of target faces and iterate over them
    target_face_pairs = itertools.combinations(induced_edge_crossings.keys(), 2)
    for face_pair in target_face_pairs:
        face_a, face_b = face_pair

        pair_induced_crossings.update({subface: dict() for subface in induced_edge_crossings[face_a].keys()})
        neighbor_assignment.update({subface: dict() for subface in induced_edge_crossings[face_a].keys()})

        # Iterate over each target face's subfaces and extract their induced edge crossings numbers
        for subface_a in induced_edge_crossings[face_a].keys():
            crossings_a = induced_edge_crossings[face_a][subface_a]

            for subface_b in induced_edge_crossings[face_b].keys():
                crossings_b = induced_edge_crossings[face_b][subface_b]

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


def select_vertex_splits(pair_induced_crossings):
    minimum_induced_crossings = float("inf")
    number_of_ties = 0
    target_subfaces = ()
    considerations = 0
    for subface_a in pair_induced_crossings.keys():
        for subface_b in pair_induced_crossings[subface_a].keys():
            considerations += 1

            if pair_induced_crossings[subface_a][subface_b] < minimum_induced_crossings:
                minimum_induced_crossings = pair_induced_crossings[subface_a][subface_b]
                target_subfaces = (subface_a, subface_b)
                number_of_ties = 0
            elif pair_induced_crossings[subface_a][subface_b] == minimum_induced_crossings:
                number_of_ties += 1
            else:
                continue

    #print(f"Considered {considerations} possible subface combinations.")
    return minimum_induced_crossings, target_subfaces, number_of_ties


def place_split_vertices(faces, selected_sub_faces, centroids, target_vertex, graph, positions):

    # Create new graph and position objects
    split_graph, split_positions = copy.deepcopy(graph), copy.deepcopy(positions)

    # Determine the indices of the to-be-split vertices
    # TODO: fix the indexing to make sense
    second_index = (max(graph.nodes()) if max(graph.nodes()) > target_vertex else target_vertex) + 1
    new_indices = [target_vertex, second_index]

    # Add newly split vertices to the graph and positions
    [split_graph.add_node(node_for_adding=split_vertex, split=1, real=1) for split_vertex in new_indices]
    sub_faces = list(selected_sub_faces.keys())
    split_positions[new_indices[0]] = centroids[faces[0]].get(sub_faces[0])
    split_positions[new_indices[1]] = centroids[faces[1]].get(sub_faces[1])

    #
    for ind, s_face in enumerate(sub_faces):

        vertices = list(selected_sub_faces[s_face])
        [split_graph.add_edge(u_of_edge=new_indices[ind], v_of_edge=vertex, real=1) for vertex in vertices]

    return split_graph, split_positions


def get_neighbors(target_face, target_subface, other_subface, neighbor_assignment, remaining_targets, all_targets):

    # Get all in-face neighbors
    neighbors = [vertex for vertex in list(target_face) if vertex in all_targets]

    # Determine which of the out-of-face neighbors belong to the current target subface
    if target_subface in neighbor_assignment.keys():
        assignments = neighbor_assignment[target_subface][other_subface]
    else:
        assignments = neighbor_assignment[other_subface][target_subface]

    # Iterate over all remaining neighbors
    for index, assignment in enumerate(assignments):
        if target_subface == assignment:
            # TODO: this set -> list conversion is unsafe as the order may not have been preserved
            neighbors.append(list(remaining_targets)[index])\

    # Return Neighbors
    return neighbors
