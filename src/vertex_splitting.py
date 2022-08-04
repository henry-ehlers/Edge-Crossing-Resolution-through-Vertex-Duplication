def planarize_graph(graph):
    #  return some graph
    return None


def locate_edge_crossings(graph, positions):
    for edge in graph.edges:
        print(edge)
        print(edge[0])
        print(edge[1])
    #  return two dicts, one for vertices and one for edges
    return None


def line_intersection(p1, p2, p3, p4):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
    if denom == 0: # parallel
        return None
    ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
    if ua < 0 or ua > 1:
        return None
    ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
    if ub < 0 or ub > 1:
        return None
    x = x1 + ua * (x2-x1)
    y = y1 + ua * (y2-y1)
    return (x, y)


# From: https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])


def sign(p1, p2, p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]) > 0


def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


def check_pairwise_edge_crossing(edge_a, edge_b, vertices):
    # Check if any two terminal vertices are identical. if so, cannot be intersecting
    if edge_a[0] == edge_b[0] or edge_a[1] == edge_b[1] or edge_a[1] == edge_b[0] or edge_a[0] == edge_b[1]:
        return False
    else:
        return intersect(vertices[edge_a[0]], vertices[edge_a[1]], vertices[edge_b[0]], vertices[edge_b[1]])


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