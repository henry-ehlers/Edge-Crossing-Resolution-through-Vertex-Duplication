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


def project_point_onto_line(point, start_point, end_point):
    normalized = (end_point-start_point) / np.linalg.norm(x=end_point-start_point)
    return np.dot(point - start_point, normalized)


def sort_vertices_along_edge(edge, vertex_set, positions):

    # Extract the 2D coordinates of the edge line
    start_vertex, end_vertex = edge[0], edge[1]
    start_point, end_point = positions[start_vertex], positions[end_vertex]

    # Initialize container for magnitudes from edge start-point
    projections = np.empty(len(vertex_set))

    # Iterate over all vertices, including start and end points of the edge
    # vertices_to_sort = list(vertex_set) + [start_vertex, end_vertex]
    for index, vertex in enumerate(vertex_set):
        projections[index] = project_point_onto_line(positions[vertex], start_point, end_point)

    # Sort indices of projections and sort vertex indices
    sorted_indices = np.argsort(projections)
    sorted_vertices = [vertex_set[i] for i in sorted_indices]

    # Ensure Start and End-Points sorted correctly (i.e first and last)
    # TODO: NUMERICAL ERROR FOR LARGE GRAPHS (WITH SMALL DIFFERENCES BETWEEN POINTS)
    # TODO: EITHER FORCE POSITIONS OR MERGE VERY CLOSELY LOCATED NODES (SEGMENT-INDUCED ONES)
    # TODO: BUILD A WHILE LOOP WHICH ITERATIVELY MERGES THE TWO CLOSEST VERTICES (PRESERVING TERMINI) UNTIL ORDER IS OK
    # assert sorted_indices[0] == (len(vertices_to_sort)-2) and sorted_indices[-1] == (len(vertices_to_sort)-1), \
    #     "Start and End Points of vector did not sort as expected"
    sorted_vertices = [start_vertex] + sorted_vertices + [end_vertex]

    # Return sorted vertex indices
    return sorted_vertices


def remove_target_vertex(graph, positions, target_vertex):

    # Remove Embedded Location of Target from position dictionary
    remaining_positions = copy.deepcopy(positions)
    del remaining_positions[target_vertex]

    # Copy Graph and remove target node from it
    remaining_graph = copy.deepcopy(graph)
    remaining_graph.remove_node(target_vertex)

    # Return copies of graph and positions without target
    return remaining_graph, remaining_positions


def get_remaining_edge_crossings(graph, edge_crossings, target_vertex):
    remaining_edge_crossings = copy.deepcopy(edge_crossings)
    edges_to_be_removed = graph.edges(nbunch=target_vertex)

    edges_removed = 0
    for edge_a in edge_crossings.keys():
        if edge_a in edges_to_be_removed:
            edges_removed += len(edge_crossings[edge_a])
            del remaining_edge_crossings[edge_a]
            continue
        for edge_b in edge_crossings[edge_a].keys():
            if edge_b in edges_to_be_removed:
                edges_removed += 1
                del remaining_edge_crossings[edge_a][edge_b]

    return remaining_edge_crossings


def remove_edges(graph, edges_to_be_removed):
    for edge in edges_to_be_removed:
        graph.remove_edge(u=edge[0], v=edge[1])


def add_virtual_edges(graph, positions, edge_to_virtual_vertex):

    # A Map of virtual edges which describe the same edge
    virtual_edge_set = {edge: [] for edge in edge_to_virtual_vertex}

    # Iterate over all edges in the graph
    for edge in edge_to_virtual_vertex.keys():

        edge_data = graph.get_edge_data(v=edge[0], u=edge[1], default={})
        # Skip edge if it does not have any edge crossings
        if len(edge_to_virtual_vertex[edge]) == 0:
            continue

        # Extract all the virtual vertices and (together with real edge points) sort them
        virtual_vertices = list(edge_to_virtual_vertex[edge])
        sorted_vertex_targets = sort_vertices_along_edge(edge, virtual_vertices, positions)

        # Connect vertices in sort order (undirected edges so order doesn't matter)
        for index in range(1, len(sorted_vertex_targets)):
            vertex_a, vertex_b = sorted_vertex_targets[index-1], sorted_vertex_targets[index]
            virtual_edge_set[edge].append((vertex_a, vertex_b))
            graph.add_edge(u_of_edge=vertex_a,
                           v_of_edge=vertex_b,
                           virtual=1,
                           real=edge_data.get("real", 0))
    return virtual_edge_set


def planarize_graph(graph, positions, edge_crossings, largest_index=None):

    # Extract basic properties of graph
    # TODO: THIS INDEX BREAK EVERYTHING SOMEHOW
    index = largest_index if largest_index is not None else max(graph.nodes)
    print(f"starting index: {index}")
    edges = list(graph.edges)  # create list for easier indexing

    edge_to_virtual_vertex = {edge: set() for edge in edges}  # have to ensure
    edges_to_be_removed = set()  # could be initialized using size of dictionary 'edge_crossings'

    # Iterate over all found edge crossings
    for edge_a in edge_crossings.keys():
        for edge_b in edge_crossings[edge_a].keys():

            # Update index
            index += 1

            # Add new vertex to graph and drawing's locations
            graph.add_node(node_for_adding=index, virtual=1)
            positions[index] = np.asarray(edge_crossings[edge_a][edge_b])

            # Log connections to new virtual vertex to be added and original (real) edges to be removed
            [edge_to_virtual_vertex[edge].add(index) for edge in [edge_a, edge_b]]
            [edges_to_be_removed.add(edge) for edge in [edge_a, edge_b]]

    # Remove original edge set and add virtual edge set
    virtual_edge_set = add_virtual_edges(graph, positions, edge_to_virtual_vertex)
    remove_edges(graph, list(edges_to_be_removed))

    #  return some new graph and new vertex positions
    return virtual_edge_set


def locate_edge_crossings(graph, positions):

    # Create object of edges for easier use
    edges = list(graph.edges)

    # Initialize vector and edge crossing containers
    vertex_crossings = {vertex: 0 for vertex in graph.nodes()}
    edge_crossings = dict()

    # ALl to all comparison of edges
    for edge_index_a in range(0, len(edges)):
        for edge_index_b in range(edge_index_a + 1, len(edges)):

            # Extract edges from edge list
            edge_a = edges[edge_index_a]
            edge_b = edges[edge_index_b]

            # Check if the two edges share a common vertex (causes numerical issues)
            if (edge_a[0] in edge_b) or (edge_a[1] in edge_b):
                continue

            # Check whether edges intersect and (if so) where
            intersection = edge_intersection(edge_a, edge_b, positions)
            if intersection is None:
                continue

            # Append edge crossing position for edges
            if edge_a not in edge_crossings:
                edge_crossings[edge_a] = dict()
            edge_crossings[edge_a][edge_b] = intersection

            # Increment edge crossing count for all vertices involves in crossing
            crossing_vertices = np.append(np.asarray(edge_a), np.asarray(edge_b))
            for vertex in crossing_vertices:
                vertex_crossings[vertex] += 1

    #  return two dicts, one for vertices and one for edge
    return edge_crossings, vertex_crossings


def is_without_edge_crossings(graph, positions):
    edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
    assert vertex_edge_crossing_equality(vertex_crossings, edge_crossings)
    return sum(vertex_crossings) == 0


def squared_distance(point_a, point_b):
    (x1, y1) = point_a
    (x2, y2) = point_b
    return (x2 - x1) ** 2 + (y2 - y1) ** 2


def find_closest_edge_intersection(edge_points, other_edges, graph, positions, must_be_real=False):
    intersections, distances = dict(), dict()
    point_a, point_b = edge_points

    for edge in other_edges:

        # Check whether the intersection is with a real edge
        if must_be_real:
            fields = graph.get_edge_data(v=edge[0], u=edge[1], default=None)
            if fields:
                if fields.get("real", 0) == 0:
                    continue
            else:
                print(f"fields of edge {edge[0], edge[1]}: {fields}")  # todo: certain edges don't have data, which shouldn't happen

        # Find
        point_c, point_d = positions[edge[0]], positions[edge[1]]
        intersection = line_intersection(point_a, point_b, point_c, point_d)
        if intersection is None: continue
        intersections[edge] = intersection
        distances[edge] = squared_distance(point_a, intersection)
    closest_intersection = min(distances, key=distances.get)

    # Return the Edge Name, and it's intersection as a tuple
    return closest_intersection, intersections[closest_intersection]


def edge_intersection(edge_a, edge_b, vertex_positions):
    point_a_0 = vertex_positions[edge_a[0]]
    point_a_1 = vertex_positions[edge_a[1]]
    point_b_0 = vertex_positions[edge_b[0]]
    point_b_1 = vertex_positions[edge_b[1]]
    return line_intersection(point_a_0, point_a_1, point_b_0, point_b_1)


def line_intersection(p1, p2, p3, p4):
    x1, y1 = float(p1[0]), float(p1[1])
    x2, y2 = float(p2[0]), float(p2[1])
    x3, y3 = float(p3[0]), float(p3[1])
    x4, y4 = float(p4[0]), float(p4[1])

    denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    if denominator == 0:  # parallel
        return None
    ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator

    # TODO: investigate these statements. just adding >= instead of > strikes me as dangerous
    if ua <= 0 or ua >= 1:
        return None
    ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator

    # TODO: investigate these statements. just adding >= instead of > strikes me as dangerous
    if ub <= 0 or ub >= 1:
        return None
    x = x1 + ua * (x2 - x1)
    y = y1 + ua * (y2 - y1)
    return x, y


def get_target_vertex(vertex_crossings, graph):
    # TODO: maybe also account for the number of edge crossings a vertex has already undergone?

    # Find all vertices which are involved in the maximal number of edge crossings
    max_score = max(vertex_crossings.values())
    potential_targets = [k for k in vertex_crossings if vertex_crossings[k] == max_score]

    # If only one vertex fulfills the maximum criterion, return said vertex's index
    if len(potential_targets) == 1:
        target_vertex = potential_targets[0]

    # If multiple vertices fulfill the criterion, return the one with the lowest degree
    else:
        adjacency = [len(graph[vertex]) for vertex in potential_targets]
        target_vertex = potential_targets[np.argmin(adjacency)]

    # Return the target vertex's index
    return target_vertex


def vertex_edge_crossing_equality(vertex_crossings, edge_crossings):
    # Calculate all vertices' total edge crossing number
    vertex_sum = 0
    for vertex in vertex_crossings.keys():
        vertex_sum += vertex_crossings[vertex]

    # Calculate all recording unique edge crossings (for non None entries)
    edge_sum = 0
    for edge_index in edge_crossings.keys():
        edge_sum += len(edge_crossings[edge_index])

    return edge_sum == (vertex_sum / 4)
