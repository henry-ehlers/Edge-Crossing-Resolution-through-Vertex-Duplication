from collections import defaultdict
import numpy as np
import copy

def sort_vertices(vertices):

    return None


def planarize_graph(graph, positions, edge_crossings):

    # Extract basic properties of graph
    n_vertices = graph.number_of_nodes()
    edges = list(graph.edges)  # create list for easier indexing

    # Initialize new, planar graph
    planar_graph = copy.deepcopy(graph)
    planar_positions = copy.deepcopy(positions)
    edge_to_virtual_vertex = {edge: set() for edge in edges}  # have to ensure
    edges_to_be_removed = set()  # could be initialized using size of dictionary 'edge_crossings'

    # Iterate over all found edge crossings
    for edge_index_a in edge_crossings.keys():
        edge_a = edges[edge_index_a]
        for edge_index_b in edge_crossings[edge_index_a].keys():
            edge_b = edges[edge_index_b]

            # Add new vertex to graph and drawing's locations
            planar_graph.add_node(node_for_adding=n_vertices, split=0, target=0, virtual=1)
            planar_positions[n_vertices] = np.asarray(edge_crossings[edge_index_a][edge_index_b])

            # Log connections to new virtual vertex to be added and original (real) edges to be removed
            [edge_to_virtual_vertex[edge].add(n_vertices) for edge in [edge_a, edge_b]]
            [edges_to_be_removed.add(edge) for edge in [edge_a, edge_b]]

            # Update index
            n_vertices += 1

    # Remove Original Edges that were involved in one or more crossing
    for edge in edges_to_be_removed:
        planar_graph.remove_edge(u=edge[0], v=edge[1])

    #  return some new graph and new vertex positions
    return planar_graph, planar_positions


def debug_edge_crossings(graph, edge_crossings):

    # Extract Edges for easier indexing
    edges = list(graph.edges)

    # ALl to all comparison of edges
    for edge_index_a in edge_crossings.keys():
        for edge_index_b in edge_crossings[edge_index_a].keys():

            # Extract edges from edge list
            edge_a = edges[edge_index_a]
            edge_b = edges[edge_index_b]

            # Print Crossings:
            print("{} and {} - {}".format(edge_a, edge_b, edge_crossings[edge_index_a][edge_index_b]))


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
            if (edge_a[0] in edge_b) or (edge_a[1] in edge_b): continue

            # Check whether edges intersect and (if so) where
            intersection = edge_intersection(edge_a, edge_b, positions)
            if intersection is None: continue

            # Append edge crossing position for edges
            if edge_index_a not in edge_crossings:
                edge_crossings[edge_index_a] = dict()
            edge_crossings[edge_index_a][edge_index_b] = intersection

            # Increment edge crossing count for all vertices involves in crossing
            crossing_vertices = np.append(np.asarray(edge_a), np.asarray(edge_b))
            for vertex_index in crossing_vertices: vertex_crossings[vertex_index] += 1

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
    if ua < 0 or ua > 1:
        return None
    ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator
    if ub < 0 or ub > 1:
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
