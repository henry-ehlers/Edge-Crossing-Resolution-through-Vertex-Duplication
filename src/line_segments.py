import numpy as np


def unbounded_line_intersection(p1, p2, p3, p4):

    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    if denominator == 0:  # parallel
        return None
    ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator

    x = x1 + ua * (x2 - x1)
    y = y1 + ua * (y2 - y1)

    return (x, y), ua


def extend_line(point_1, point_2, bounds=((-2, -2), (-2, 2), (2, 2), (2, -2))):
    print(bounds)
    ts = np.array(object=[float("-inf"), float("inf")], dtype=float)
    intersections = [(None, None), (None, None)]
    for index_a in range(0, len(bounds)):
        index_b = (index_a + 1) % len(bounds)
        result = unbounded_line_intersection(point_1, point_2, bounds[index_a], bounds[index_b])
        if result is None:
            continue
        print("t: {}".format(result[1]))
        print("intersection: {}".format(result[0]))
        if 0 >= result[1] >= ts[0]:
            ts[0] = result[1]
            intersections[0] = result[0]
        elif 1 <= result[1] <= ts[1]:
            ts[1] = result[1]
            intersections[1] = result[0]

    print(f"found intersections: {intersections}")
    return intersections

