from src.edge_crossings import *
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
    print(f"real edges: {real_edges}")

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

    # TODO: turn induced edge crossings into table
    selected_face = list(centroids.keys())[0]
    get_subface_edge_crossing_table(face=selected_face,
                                    subface_crossings=induced_edge_crossings[selected_face],
                                    target_vertices=target_neighbors)
    # Return Complete Set
    return induced_edge_crossings


def get_subface_edge_crossing_table(face, subface_crossings, target_vertices):

    # Create Lists of sub-faces and faces
    sub_faces = [subface for subface in subface_crossings.keys()]
    faces = [face] * len(subface_crossings)

    # Tabulate the number of crossings for all target vertex and subface combinations
    crossings = np.empty(shape=(len(subface_crossings), len(target_vertices)),
                         dtype=int)
    for index, subface in enumerate(sub_faces):
        print(f"subface: {subface}")
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


def place_split_vertices(graph, positions, target_vertex, all_targets, remaining_targets,
                         neighbor_assignment, centroids, target_subfaces, face_subface_map):

    # Create new graph and position objects
    split_graph, split_positions = copy.deepcopy(graph), copy.deepcopy(positions)

    # Determine the indices of the to-be-split vertices
    split_vertex_indices = [target_vertex, max(graph.nodes()) + 1]

    # Neighbors connected
    neighbors_connected = []

    # Iterate over all selected subfaces
    for subface_index, target_subface in enumerate(target_subfaces):
        #print(f"Subface #{subface_index} - {target_subface}")

        # Determine which super face the subface falls within
        target_face = [face for face in face_subface_map.keys() if target_subface in face_subface_map[face]][0]
        #print(set(target_subfaces))
        other_subface = target_subfaces[0] if target_subfaces[1] == target_subface else target_subfaces[1]
        #print(f"Target Face: {target_face}")
        #print(f"Other Subface: {other_subface}")

        # Neighborhood
        neighbors = get_neighbors(target_face, target_subface, other_subface, neighbor_assignment, remaining_targets, all_targets)
        #print(f"Neighbors: {neighbors}")

        # Extract the location of the split vertex
        location = centroids[target_face][target_subface]
        vertex = split_vertex_indices[subface_index]
        #print(f"Vertex #{vertex} at location {location}")

        # Place Split vertex
        #print(f"Adding Vertex {vertex}")
        split_graph.add_node(node_for_adding=vertex, split=1, target=0, virtual=0, boundary=0, segment=0)
        split_positions[vertex] = np.array(location)

        # Add Edges
        for neighbor_vertex in neighbors:
            split_graph.add_edge(u_of_edge=vertex, v_of_edge=neighbor_vertex, virtual=0, target=0, segment=0)
            neighbors_connected.append(neighbor_vertex)

    return split_graph, split_positions, neighbors_connected


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
