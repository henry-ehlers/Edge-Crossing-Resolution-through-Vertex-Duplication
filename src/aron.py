import numpy as np

def cross(a, b):
    return a[0] * b[1] - b[0] * a[1]

def vector_angle(vector_1, vector_2):
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)

    dot_product = np.dot(unit_vector_1, unit_vector_2)
    # print(vector_1, vector_2, unit_vector_1, unit_vector_2, dot_product)
    angle = np.arccos(dot_product)

    #angle = np.arcsin(cross(unit_vector_1, unit_vector_2))
    #if angle < 0:
    #    angle *= -1
    #else:
    #    angle = 2 * np.pi - angle

    return np.degrees(angle)


# cw = negative
# ccw = positive
def area(x, y):
    return sum(x[i] * (y[i + 1] - y[i - 1]) for i in range(-1, len(x) - 1)) / 2.0

x = [0, 1, 0]
y = [0, 0, 1]
print(area(x, y))


edges = [np.array(point) for point in [(0, 0), (1, 1), (2, 1), (3, 0), (4, 1), (5, 4), (3, 3), (2, 2), (1, 2), (0, 3), (1, 4), (-1, 4), (0, 1)]]
#edges = [np.array(point) for point in [(0, 0), (1, 1), (2, 1), (3, 0), (3, 3), (0, 3)]]
#edges = [np.array(point) for point in reversed([(0, 0), (1, 1), (2, 1), (3, 0), (3, 3), (0, 3)])]
print(edges)
for angle_around in range(0, len(edges)):
    vector_1 = edges[angle_around - 1] - edges[angle_around]
    vector_2 = edges[(angle_around + 1) % len(edges)] - edges[angle_around]
#    print(f"Angle Around {edges[angle_around]}: {vector_angle(vector_1, vector_2)}")

for i in range(len(edges)):
    p1 = edges[i]
    ref = edges[i - 1]
    p2 = edges[i - 2]
    a = p1 - ref
    b = p2 - ref

    print('Point', ref)

    if cross(a, b) > 0.0:
        print('Angle', vector_angle(a, b))
    else:
        print('Angle', 360 - vector_angle(a, b))
    print('')
