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

    print("EDGE SETS:")
    [print(edge_set) for edge_set in virtual_edge_set]
    print("---------------------------")

    # Store nodes and edges for easier look-up
    node_list = list(segment_graph.nodes())
    edges = frozenset([frozenset(edge) for edge in list(segment_graph.edges())])

    # Store the number of nodes and largest index
    vertex_index = max(node_list)
    number_of_nodes = len(node_list)

    [print(edge_set) for edge_set in virtual_edge_set]

    # Iterate over all pairwise vertex combinations
    for index_a in range(0, number_of_nodes):
        vertex_a = node_list[index_a]
        for index_b in range(index_a + 1, number_of_nodes):
            vertex_b = node_list[index_b]

            # Check if this particular combination of vertices maps to a virtual edge set
            if any([{vertex_a, vertex_b} <= v_edge_set for v_edge_set in virtual_edge_set]):
                # Ensure that both vertices are real
                if any(graph.nodes[vertex]["virtual"] == 1 for vertex in [vertex_a, vertex_b]):
                    continue
                else:
                    already_connected = 1
            else:
                # Check whether an edge already exists between vertices a and b
                if {node_list[index_a], node_list[index_b]} in edges:
                    already_connected = 1
                else:
                    already_connected = 0

            # Calculate intersections with boundary
            intersections = extend_line(segment_positions[vertex_a], segment_positions[vertex_b], bounds)

            # Add new boundary vertices to the graph
            for added_vertex in range(0, len(intersections)):
                vertex_index += 1
                segment_graph.add_node(node_for_adding=vertex_index, split=0, target=0, virtual=0, boundary=1)
                segment_positions[vertex_index] = np.asarray(intersections[added_vertex])

            # Determine edge connections based on whether an edge already exists between vertex A and B
            if already_connected:
                connections = [(vertex_index-1, vertex_a), (vertex_b, vertex_index)]
            else:
                connections = [(vertex_index-1, vertex_a), (vertex_a, vertex_b), (vertex_b, vertex_index)]

            # Add new edges (line segments) to graph
            for connection in connections:
                segment_graph.add_edge(u_of_edge=connection[0], v_of_edge=connection[1], virtual=0, target=0, segment=1)

    # Return new graph and positions objects
    return segment_graph, segment_positions


def cull_all_line_segment_graph(graph, positions, target_faces, face_edge_map):
    culled_graph = copy.deepcopy(graph)
    culled_positions = copy.deepcopy(positions)
    # target_face_edges = [face_edge_map[target_face] for target_face in target_faces]

    # Keep track of non-intersecting segments
    edges_to_be_removed = set()
    nodes_to_be_removed = set()

    # Initialize empty nested dictionary
    face_intersection_map = {target_face: dict() for target_face in target_faces}
    print(f"face intersection map: {face_intersection_map}")
    print(f"target faces: {target_faces}")

    for edge in culled_graph.edges:
        vertex_a, vertex_b = edge[0], edge[1]

        if culled_graph.edges[edge]["segment"]:
            intersection_found = False

            for target_face in target_faces:
                print(f"target face: {target_face}")
                for face_edge in face_edge_map[target_face]:
                    #print(f"face edge: {face_edge}")

                    # Define Points and calculate intersection between them
                    edge_point_a, edge_point_b = culled_positions[vertex_a], culled_positions[vertex_b]
                    face_point_a, face_point_b = culled_positions[face_edge[0]], culled_positions[face_edge[1]]
                    intersection = line_intersection(edge_point_a, edge_point_b, face_point_a, face_point_b)
                    print(f"Intersection: {intersection}")

                    # TODO: keep track of intersections and add new virtual vertices
                    # TODO: store intersections as dictionary of lists (frozenset as key)
                    # TODO: don't check edges that are in the the face edge set

                    # Store intersection if one exists
                    if intersection is not None:
                        intersection_found = True
                        if (vertex_a, vertex_b) not in face_intersection_map[target_face].keys():
                            face_intersection_map[target_face][(vertex_a, vertex_b)] = dict()
                        face_intersection_map[target_face][(vertex_a, vertex_b)][face_edge] = intersection

            if not intersection_found:
                edges_to_be_removed.add(edge)
                if any([culled_graph.nodes[vertex]["boundary"] == 1 for vertex in [vertex_a, vertex_b]]):
                    nodes_to_be_removed.add(vertex_a if culled_graph.nodes[vertex_a]["boundary"] == 1 else vertex_b)

    for edge in edges_to_be_removed:
        culled_graph.remove_edge(u=edge[0], v=edge[1])
    for node in nodes_to_be_removed:
        culled_graph.remove_node(node)

    return culled_graph, culled_positions, face_intersection_map


def create_subface_graph(graph, positions, target_faces, face_intersection_map):

    # Iterate over all (hopefully 2) target faces
    for target_face in target_faces:
        for intersecting_edge in face_intersection_map[target_face].keys():
            print(f"Intersecting Edge: {intersecting_edge}")
            if len(intersecting_edge.keys()) > 2:
                print("whoops!")

    return None
