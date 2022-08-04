import numpy as np


def planarize_graph(graph):
    #  return some graph
    return None


def locate_edge_crossings(graph, positions):

    # Initialize vector and edge crossing containers
    edges = list(graph.edges)
    vertex_crossings = np.zeros(shape=len(graph.nodes), dtype=int)

    # ALl to all comparison of edges
    for edge_index_a in range(len(edges)):
        for edge_index_b in range(edge_index_a+1, len(edges)):
            print("{} - {}".format(edge_index_a, edge_index_b))
            # Skip if the two edges are identical
            if edge_index_a == edge_index_b: continue

            # Extract edges from edge list
            edge_a = edges[edge_index_a]
            edge_b = edges[edge_index_b]

            # Check whether edges intersect and (if so) where
            intersection = edge_intersection(edge_a, edge_b, positions)
            if intersection is None: continue

            # Append edge crossing information
            for vertex_index in np.append(np.asarray(edge_a), np.asarray(edge_b)):
                vertex_crossings[vertex_index] += 1

    print(vertex_crossings)
    #  return two dicts, one for vertices and one for edge
    return None


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
    denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    if denom == 0:  # parallel
        return None
    ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denom
    if ua < 0 or ua > 1:
        return None
    ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom
    if ub < 0 or ub > 1:
        return None
    x = x1 + ua * (x2 - x1)
    y = y1 + ua * (y2 - y1)
    return (x, y)



# # Investigate Edges (upper triangle)
# edges = list(G.edges)
# intersections = list()
# for edge_a in range(len(edges)):
#     for edge_b in range(edge_a + 1, len(edges)):
#         if check_pairwise_edge_crossing(edges[edge_a], edges[edge_b], pos):
#             intersections.append([edges[edge_a], edges[edge_b]])
#
# # Calculate Crossing Number
# crossing_number = dict()
# for intersection in intersections:
#     for edge in intersection:
#         if edge in crossing_number:
#             crossing_number[edge] += 1
#         else:
#             crossing_number[edge] = 1
# print("Crossing Numbers:")
# print(crossing_number)
#
# # find edge with most crossings
# maxima = find_maxima(crossing_number)
#
# # Find points of intersection along split-target edge
# # TODO: hunting down a bug in the crossing number
# crossings = list()
# for intersection in intersections:
#     if split_target in intersection:
#         if split_target == intersection[0]:
#             crossing_edge = intersection[1]
#         else:
#             crossing_edge = intersection[0]
#         a1 = pos[split_target[0]]
#         a2 = pos[split_target[1]]
#         b1 = pos[crossing_edge[0]]
#         b2 = pos[crossing_edge[1]]
#         crossings.append(line_intersection(a1, a2, b1, b2))
