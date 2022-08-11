from src.edge_crossings import *
import numpy as np
import copy


def unbounded_line_intersection(p1, p2, p3, p4):

    # Extract individual floats from point tuples
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    # Calculate denominator and return None if lines are parallel
    denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    if denominator == 0:
        return None

    # Calculate scalar
    ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator

    # Calculate new point coordinates
    x = x1 + ua * (x2 - x1)
    y = y1 + ua * (y2 - y1)

    # Return Coordinates and scalar
    return (x, y), ua


def extend_line(point_1, point_2, bounds=((-2, -2), (-2, 2), (2, 2), (2, -2))):

    ts = np.array(object=[float("-inf"), float("inf")], dtype=float)
    intersections = [(None, None), (None, None)]
    for index_a in range(0, len(bounds)):
        index_b = (index_a + 1) % len(bounds)
        result = unbounded_line_intersection(point_1, point_2, bounds[index_a], bounds[index_b])
        if result is None:
            continue
        if 0 >= result[1] >= ts[0]:
            ts[0] = result[1]
            intersections[0] = result[0]
        elif 1 <= result[1] <= ts[1]:
            ts[1] = result[1]
            intersections[1] = result[0]

    return intersections


def draw_all_line_segments(graph, positions, virtual_edge_set, bounds=((-1, -1), (-1, 1), (1, 1), (1, -1))):

    # Create new graph and positions objects
    segment_graph = copy.deepcopy(graph)
    segment_positions = copy.deepcopy(positions)

    new_virtual_edge_sets = set()

    # Store nodes and edges for easier look-up
    node_list = list(segment_graph.nodes())
    edges = frozenset([frozenset(edge) for edge in list(segment_graph.edges())])

    # Store the number of nodes and largest index
    vertex_index = max(node_list)
    number_of_nodes = len(node_list)

    # Iterate over all pairwise vertex combinations
    for index_a in range(0, number_of_nodes):
        vertex_a = node_list[index_a]
        for index_b in range(index_a + 1, number_of_nodes):
            vertex_b = node_list[index_b]

            # Check if this particular combination of vertices maps to a virtual edge set
            found_set = [v_edge_set for v_edge_set in virtual_edge_set if {vertex_a, vertex_b} <= v_edge_set]
            if found_set:
                if any(graph.nodes[vertex]["virtual"] == 1 for vertex in [vertex_a, vertex_b]):
                    continue
                already_connected = 1
            else:
                already_connected = 1 if {node_list[index_a], node_list[index_b]} in edges else 0

            # Calculate intersections with boundary
            intersections = extend_line(segment_positions[vertex_a], segment_positions[vertex_b], bounds)

            # Add new boundary vertices to the graph
            for added_vertex in range(0, len(intersections)):
                vertex_index += 1
                segment_graph.add_node(node_for_adding=vertex_index, split=0, target=0, virtual=0, boundary=1, segment=0)
                segment_positions[vertex_index] = np.asarray(intersections[added_vertex])

            if found_set:
                new_set = frozenset(list(found_set[0])+[vertex_index - 1, vertex_index])
                new_virtual_edge_sets.add(frozenset(new_set))

            # Add new edges (line segments) to graph
            connections = [(vertex_index - 1, vertex_a), (vertex_b, vertex_index)] if already_connected \
                else [(vertex_index - 1, vertex_a), (vertex_a, vertex_b), (vertex_b, vertex_index)]
            for connection in connections:
                segment_graph.add_edge(u_of_edge=connection[0], v_of_edge=connection[1], virtual=0, target=0, segment=1)

    # Return new graph and positions objects
    return segment_graph, segment_positions, new_virtual_edge_sets


def cull_all_line_segment_graph(graph, positions, target_faces, face_edge_map, virtual_edge_set):

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

    for target_face in target_faces:
        for intersecting_edge in face_intersection_map[target_face].keys():
            edges_to_be_removed.add(frozenset(intersecting_edge))
            edge_targets = []
            for face_edge in face_intersection_map[target_face][intersecting_edge].keys():
                if type(face_edge) is tuple:
                    if face_edge not in edge_to_virtual_vertex:
                        edge_to_virtual_vertex[face_edge] = set()
                    edges_to_be_removed.add(frozenset(face_edge))
                    intersection = face_intersection_map[target_face][intersecting_edge][face_edge]
                    vertex_index += 1
                    edge_targets.append(vertex_index)
                    graph.add_node(node_for_adding=vertex_index, split=0, target=0, virtual=0, boundary=0, segment=1)
                    positions[vertex_index] = np.asarray(intersection)
                    edge_to_virtual_vertex[face_edge].add(vertex_index)
                else:
                    edge_targets.append(face_edge)
            graph.add_edge(u_of_edge=edge_targets[0], v_of_edge=edge_targets[1], virtual=0, target=0, segment=1)

    # Add virtual edge connections
    [print(f"{index} - {edge_to_virtual_vertex[index]}") for index in edge_to_virtual_vertex.keys()]
    virtual_edge_set = add_virtual_edges(graph, positions, edge_to_virtual_vertex)

    # Remove
    [graph.remove_edge(u=list(edge)[0], v=list(edge)[1]) for edge in edges_to_be_removed]
    [graph.remove_node(vertex) for vertex in nodes_to_be_removed]

    return graph, positions, virtual_edge_set

