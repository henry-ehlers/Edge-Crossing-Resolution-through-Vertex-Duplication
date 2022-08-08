import numpy as np
import copy


def project_point_onto_line(point, start_point, end_point):
    normalized = (end_point-start_point) / np.linalg.norm(x=end_point-start_point)
    return np.dot(point - start_point, normalized)


def sort_vertices_along_edge(edge, vertex_set, positions):

    # Extract the 2D coordinates of the edge line
    start_vertex, end_vertex = edge[0], edge[1]
    start_point, end_point = positions[start_vertex], positions[end_vertex]

    # Initialize container for magnitudes from edge start-point
    projections = np.empty(len(vertex_set)+2)

    # Iterate over all vertices, including start and end points of the edge
    vertices_to_sort = list(vertex_set) + [start_vertex, end_vertex]
    print("vertices_to_sort: {}".format(vertices_to_sort))
    for index, vertex in enumerate(vertices_to_sort):
        projections[index] = project_point_onto_line(positions[vertex], start_point, end_point)
    print(projections)

    # Sort indices of projections and sort vertex indices
    sorted_indices = np.argsort(projections)
    sorted_vertices = [vertices_to_sort[i] for i in sorted_indices]
    print("sorted: {}".format(sorted_vertices))
    # Ensure Start and End-Points sorted correctly (i.e first and last)
    assert sorted_indices[0] == (len(vertices_to_sort)-2) and sorted_indices[-1] == (len(vertices_to_sort)-1), \
        "Start and End Points of vector did not sort as expected"

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
    print(edges_to_be_removed)
    print("To be deleted: {}".format(len(edges_to_be_removed)))
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

    print("Crossings Removed: {}".format(edges_removed))
    return remaining_edge_crossings


def remove_edges(graph, edges_to_be_removed):
    for edge in edges_to_be_removed:
        graph.remove_edge(u=edge[0], v=edge[1])


def add_virtual_edges(graph, positions, edge_to_virtual_vertex):

    # Iterate over all edges in the graph
    for edge in edge_to_virtual_vertex.keys():
        print("\nEdge: {}".format(edge))

        # Skip edge if it does not have any edge crossings
        if len(edge_to_virtual_vertex[edge]) == 0:
            continue

        # Extract all the virtual vertices and (together with real edge points) sort them
        virtual_vertices = list(edge_to_virtual_vertex[edge])
        sorted_vertex_targets = sort_vertices_along_edge(edge, virtual_vertices, positions)

        # Connect vertices in sort order (undirected edges so order doesn't matter)
        for index in range(1, len(sorted_vertex_targets)):
            vertex_a, vertex_b = sorted_vertex_targets[index-1], sorted_vertex_targets[index]
            graph.add_edge(u_of_edge=vertex_a, v_of_edge=vertex_b, virtual=1)


def planarize_graph(graph, positions, edge_crossings, starting_index):

    # Extract basic properties of graph
    index = starting_index
    edges = list(graph.edges)  # create list for easier indexing

    # Initialize new, planar graph
    planar_graph = copy.deepcopy(graph)
    planar_positions = copy.deepcopy(positions)
    print("Length Planar Positions: {}".format(len(planar_positions)))

    edge_to_virtual_vertex = {edge: set() for edge in edges}  # have to ensure
    edges_to_be_removed = set()  # could be initialized using size of dictionary 'edge_crossings'

    # Iterate over all found edge crossings
    for edge_a in edge_crossings.keys():
        for edge_b in edge_crossings[edge_a].keys():
            print("{} - {} : {}".format(edge_a, edge_b, index))

            # Add new vertex to graph and drawing's locations
            planar_graph.add_node(node_for_adding=index, split=0, target=0, virtual=1)
            planar_positions[index] = np.asarray(edge_crossings[edge_a][edge_b])

            # Log connections to new virtual vertex to be added and original (real) edges to be removed
            [edge_to_virtual_vertex[edge].add(index) for edge in [edge_a, edge_b]]
            [edges_to_be_removed.add(edge) for edge in [edge_a, edge_b]]

            # Update index
            index += 1

    # Remove original edge set and add virtual edge set
    add_virtual_edges(planar_graph, planar_positions, edge_to_virtual_vertex)
    remove_edges(planar_graph, list(edges_to_be_removed))

    #  return some new graph and new vertex positions
    return planar_graph, planar_positions


def locate_edge_crossings(graph, positions):

    # Create object of edges for easier use
    edges = list(graph.edges)

    # Initialize vector and edge crossing containers
    vertex_crossings = np.zeros(shape=len(graph.nodes), dtype=int)
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
            for vertex_index in crossing_vertices:
                vertex_crossings[vertex_index] += 1

    #  return two dicts, one for vertices and one for edge
    return edge_crossings, vertex_crossings


def is_without_edge_crossings(graph, positions):
    edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
    assert vertex_edge_crossing_equality(vertex_crossings, edge_crossings)
    return sum(vertex_crossings) == 0


def edge_intersection(edge_a, edge_b, vertex_positions):
    point_a_0 = vertex_positions[edge_a[0]]
    point_a_1 = vertex_positions[edge_a[1]]
    point_b_0 = vertex_positions[edge_b[0]]
    point_b_1 = vertex_positions[edge_b[1]]
    return line_intersection(point_a_0, point_a_1, point_b_0, point_b_1)


def line_intersection(p1, p2, p3, p4):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    if denominator == 0:  # parallel
        return None
    ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator

    # TODO: investigate these statements. just adding >= instead of > strikes me as dangerous
    if ua < 0 or ua >= 1:
        return None
    ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator

    # TODO: investigate these statements. just adding >= instead of > strikes me as dangerous
    if ub < 0 or ub >= 1:
        return None
    x = x1 + ua * (x2 - x1)
    y = y1 + ua * (y2 - y1)
    return x, y


def get_target_vertex_index(vertex_crossings, graph):
    # TODO: maybe also account for the number of edge crossings a vertex has already undergone?

    # Find all vertices which are involved in the maximal number of edge crossings
    potential_targets = np.argwhere(vertex_crossings == np.amax(vertex_crossings)).flatten().tolist()

    # If only one vertex fulfills the maximum criterion, return said vertex's index
    if len(potential_targets) == 1:
        target = potential_targets[0]

    # If multiple vertices fulfill the criterion, return the one with the lowest degree
    else:
        adjacency = [len(graph[target_index]) for target_index in potential_targets]
        target = potential_targets[np.argmin(adjacency)]

    # Return the target vertex's index
    return target


def vertex_edge_crossing_equality(vertex_crossings, edge_crossings):
    # Calculate all vertices' total edge crossing number
    vertex_sum = sum(vertex_crossings)

    # Calculate all recording unique edge crossings (for non None entries)
    edge_sum = 0
    for edge_index in edge_crossings.keys():
        edge_sum += len(edge_crossings[edge_index])

    return edge_sum == (vertex_sum / 4)
