import numpy as np


def planarize_graph(graph, positions, edge_crossings):

    # TODO: remove edges involved in edge crossing
    # TODO: add new vertices to drawing where edge crossings are located
    # TODO: add new edges between old and new vertices

    #  return some new graph and new vertex positions
    return None


def locate_edge_crossings(graph, positions):

    # Create object of edges for easier use
    edges = list(graph.edges)

    # Initialize vector and edge crossing containers
    vertex_crossings = np.zeros(shape=len(graph.nodes), dtype=int)
    edge_crossings = np.empty(shape=len(edges), dtype=object)

    # ALl to all comparison of edges
    for edge_index_a in range(0, len(edges)):
        for edge_index_b in range(edge_index_a+1, len(edges)):

            # Extract edges from edge list
            edge_a = edges[edge_index_a]
            edge_b = edges[edge_index_b]

            # Check if the two edges share a common vertex (causes numerical issues)
            if (edge_a[0] in edge_b) or (edge_a[1] in edge_b): continue

            # Check whether edges intersect and (if so) where
            intersection = edge_intersection(edge_a, edge_b, positions)
            if intersection is None: continue

            # Append edge crossing position for edges
            if edge_crossings[edge_index_a] is None:
                edge_crossings[edge_index_a] = [intersection]
            else:
                edge_crossings[edge_index_a].append(intersection)

            # Increment edge crossing count for all vertices involves in crossing
            crossing_vertices = np.append(np.asarray(edge_a), np.asarray(edge_b))
            for vertex_index in crossing_vertices: vertex_crossings[vertex_index] += 1

    #  return two dicts, one for vertices and one for edge
    return edge_crossings, vertex_crossings


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


def get_target_vertex_index(vertex_crossings):
    return np.argmax(vertex_crossings)


def vertex_edge_crossing_equality(vertex_crossings, edge_crossings):

    # Calculate all vertices' total edge crossing number
    vertex_sum = sum(vertex_crossings)

    # Calculate all recording unique edge crossings (for non None entries)
    edge_sum = 0
    for edge_crossing in edge_crossings:
        if not edge_crossing is None:
            edge_sum += len(edge_crossing)

    return edge_sum == (vertex_sum / 4)

