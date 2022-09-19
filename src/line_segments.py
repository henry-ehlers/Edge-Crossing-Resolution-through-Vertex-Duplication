from src.edge_crossings import *
from src.edges import *

import networkx as nx
import numpy as np
import copy


def get_nested_key(dictionary, key_a, key_b):
    try:
        return dictionary[key_a][key_b]
    except KeyError:
        return None


def draw_all_line_segments(graph, positions, virtual_edge_set, bounds, already_extended, start_index=None):

    # Create new graph and positions objects
    segment_graph, segment_positions = copy.deepcopy(graph), copy.deepcopy(positions)

    # Store nodes and edges for easier look-up
    virtual_nodes = nx.get_node_attributes(graph, "virtual")
    corner_nodes = nx.get_node_attributes(graph, "corner")

    real_nodes = [v for v in segment_graph.nodes if (virtual_nodes.get(v, 0) != 1) and (corner_nodes.get(v, 0) != 1)]
    edges = frozenset([frozenset(edge) for edge in list(segment_graph.edges())])

    # Store the number of nodes and largest index
    vertex_index = start_index if start_index else max(graph.nodes)
    number_of_nodes = len(real_nodes)

    # Iterate over all pairwise vertex combinations
    for index_a in range(0, number_of_nodes):
        vertex_a = real_nodes[index_a]
        print(f"\nVertex A: {vertex_a}")
        for index_b in range(index_a + 1, number_of_nodes):
            vertex_b = real_nodes[index_b]
            print(f"Vertex B: {vertex_b}")

            # Calculate intersections with boundary
            intersections = extend_line(segment_positions[vertex_a], segment_positions[vertex_b], bounds)

            # # Check if this particular combination of vertices is connected (by virtual edge sets)
            virtual_connection = [v_edge_set for v_edge_set in virtual_edge_set if {vertex_a, vertex_b} <= v_edge_set]
            if virtual_connection:
                already_connected = 1
            else:
                already_connected = 1 if {real_nodes[index_a], real_nodes[index_b]} in edges else 0
            print(f"already connected: {already_connected}")

            if already_connected == 0:
                segment_graph.add_edge(u_of_edge=vertex_a, v_of_edge=vertex_b, segment=1)

            for pair_index, (target, joint) in enumerate([(vertex_b, vertex_a), (vertex_a, vertex_b)]):
                print(f"target: {target}")
                print(f"joint:  {joint}")
                print(f"intersection: {intersections[pair_index]}")

                extended_vertex = get_nested_key(already_extended, target, joint)
                print(f"already extended: {extended_vertex}")

                # TODO: check if edges still/already exist

                # Check if vertex pair has already been extended
                if extended_vertex is not None:
                    extended_edges = virtual_edge_set.get(frozenset((joint, extended_vertex[0])), None)

                    # Check if any virtual edges exist between joint and extension terminus
                    if extended_edges is not None:
                        connections = get_vertex_sequence(edges=extended_edges, first_node=joint, is_ordered=False)
                        print(f"virtual edge set: {extended_edges}")
                        print(f"vertex chain: {connections}")

                    # If not, then the connection is a direct one
                    else:
                        connections = [joint, extended_vertex[0]]
                        print(f"Direct connection between {joint} and {extended_vertex[0]}")

                # If not, then create new vertex
                else:
                    connections = [joint]

                # Check if extended vertex is a boundary -> if not, create a new one
                if (extended_vertex is None) or (extended_vertex[1] != 1):
                    vertex_index += 1
                    segment_graph.add_node(node_for_adding=vertex_index, boundary=1)
                    segment_positions[vertex_index] = np.asarray(intersections[pair_index])
                    connections.append(vertex_index)
                    print(f"Added vertex {vertex_index}")

                print(f"connections: {connections}")
                # Add new edges (line segments) to graph
                for index in range(1, len(connections)):
                    print(f"index: {index}")
                    edge = {connections[index - 1], connections[index]}
                    print(f"Edge: {edge}")
                    if edge in edges:
                        continue
                    segment_graph.add_edge(u_of_edge=int(edge.pop()), v_of_edge=int(edge.pop()), segment=1)

    # Return new graph and positions objects
    return segment_graph, segment_positions


def cull_all_line_segment_graph(graph, positions, target_faces, face_edge_map):

    # Create new graph objects
    culled_graph = copy.deepcopy(graph)
    culled_positions = copy.deepcopy(positions)

    # Keep track of non-intersecting segments
    edges_to_be_removed = set()
    nodes_to_be_removed = set()

    # Initialize empty nested dictionary
    face_intersection_map = {target_face: dict() for target_face in target_faces}

    for edge in culled_graph.edges:
        vertex_a, vertex_b = edge[0], edge[1]

        if not culled_graph.edges[edge]["segment"]:
            continue

        # Keep Track of intersections found
        delete_edge = True

        for target_face in target_faces:
            intersections_found = 0

            for face_edge in face_edge_map[target_face]:

                # Determine if line actually passes through one of the vertices which define the face
                if {vertex_a, vertex_b} & {face_edge[0], face_edge[1]}:
                    intersection_key = face_edge[0] if face_edge[0] in {vertex_a, vertex_b} else face_edge[1]
                    intersection = culled_positions[intersection_key]
                else:
                    intersection_key = (face_edge[0], face_edge[1])
                    edge_point_a, edge_point_b = culled_positions[vertex_a], culled_positions[vertex_b]
                    face_point_a, face_point_b = culled_positions[face_edge[0]], culled_positions[face_edge[1]]
                    intersection = line_intersection(edge_point_a, edge_point_b, face_point_a, face_point_b)

                # Store intersection if one exists
                if intersection is not None:
                    if edge not in face_intersection_map[target_face].keys():
                        face_intersection_map[target_face][edge] = dict()
                    if intersection_key not in face_intersection_map[target_face][edge].keys():
                        intersections_found += 1
                        face_intersection_map[target_face][edge][intersection_key] = intersection

            # If only one intersection was found, remove the edge
            if intersections_found == 1:
                intersections_found = 0
                face_intersection_map[target_face].pop(edge)

            # If the edge did intersect (at least) a face, unmark it for deletion
            if intersections_found > 0:
                delete_edge = False

        # If Edge does not intersect at all, remove it
        if delete_edge:
            edges_to_be_removed.add(edge)

        # If any of the nodes part of the current edge are virtual boundary node, delete them
        if delete_edge and any([culled_graph.nodes[vertex]["boundary"] == 1 for vertex in [vertex_a, vertex_b]]):
            nodes_to_be_removed.add(vertex_a if culled_graph.nodes[vertex_a]["boundary"] == 1 else vertex_b)

    # Remove Edges and Vertices which did not intersect any face
    for edge in edges_to_be_removed:
        culled_graph.remove_edge(u=edge[0], v=edge[1])
    for node in nodes_to_be_removed:
        culled_graph.remove_node(node)

    # Return new graph, positions, and intersection map
    return culled_graph, culled_positions, face_intersection_map


def create_subface_graph(graph, positions, target_faces, face_intersection_map):

    # Iterate over all (hopefully 2) target faces
    node_list = list(graph.nodes())
    vertex_index = max(node_list)

    edges_to_be_removed = set()
    nodes_to_be_removed = [vertex for vertex in graph if graph.nodes[vertex]["boundary"] == 1]
    edge_to_virtual_vertex = dict()
    face_vertex_map = {face: set() for face in target_faces}

    for target_face in target_faces:
        [face_vertex_map[target_face].add(face_vertex) for face_vertex in list(target_face)]

        for intersecting_edge in face_intersection_map[target_face].keys():
            edges_to_be_removed.add(frozenset(intersecting_edge))
            edge_targets = []
            for face_edge in face_intersection_map[target_face][intersecting_edge].keys():

                # Ensure the intersection is with a line, not a face vertex
                if type(face_edge) is tuple:

                    # Add currently considered face edge to the dictionary
                    if face_edge not in edge_to_virtual_vertex:
                        edge_to_virtual_vertex[face_edge] = set()
                    edges_to_be_removed.add(frozenset(face_edge))
                    intersection = face_intersection_map[target_face][intersecting_edge][face_edge]

                    # Create new Vertex
                    vertex_index += 1
                    edge_targets.append(vertex_index)
                    graph.add_node(node_for_adding=vertex_index, segment=1)
                    face_vertex_map[target_face].add(vertex_index)
                    positions[vertex_index] = np.asarray(intersection)
                    edge_to_virtual_vertex[face_edge].add(vertex_index)
                else:
                    edge_targets.append(face_edge)
            graph.add_edge(u_of_edge=edge_targets[0], v_of_edge=edge_targets[1], segment=1)

    # Add virtual edge connections
    # [print(f"{index} - {edge_to_virtual_vertex[index]}") for index in edge_to_virtual_vertex.keys()]
    virtual_edge_set = add_virtual_edges(graph, positions, edge_to_virtual_vertex)

    # Remove
    [graph.remove_edge(u=list(edge)[0], v=list(edge)[1]) for edge in edges_to_be_removed]
    [graph.remove_node(vertex) for vertex in nodes_to_be_removed]

    return graph, positions, virtual_edge_set, face_vertex_map

