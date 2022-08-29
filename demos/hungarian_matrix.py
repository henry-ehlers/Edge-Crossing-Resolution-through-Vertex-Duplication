import numpy as np

# Create an n * m matrix, of n=7 workers (i.e. sight cells) and m=4 tasks (i.e. target neighbor vertices)
cost_matrx = np.array(
    [
        [0, 0, 1, 1],
        [1, 0, 0, 0],
        [0, 1, 1, 1],
        [1, 1, 0, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 0],
        [1, 0, 1, 0]
    ]
)

print(cost_matrx)
