import numpy as np

# Terminology ----------------------------------------------------------------------------------------------------------
"""
m = the number of tasks = the number of target vertex neighbors that we need to connect to one of the split vertices
In this example, m = 4, i.e. the number of columns of our two cost matrices

n_a = the number of subfaces identified in selected face (or sight cell) A. 
Here, n_a = 7, i.e the number of rows in one of our cost matrix induced_edge_crossings_a.

n_b = the number of subfaces identified in selected face (or sight cell) B. 
Here, n_a = 4, i.e the number of rows in one of our cost matrix induced_edge_crossings_a.
"""

# Some Notes -----------------------------------------------------------------------------------------------------------
"""
n_a != n_b, i.e. the number of subfaces can be different per face. 

Each cost matrix will have at least one column that is equal to zero for every row. Subfaces are created from a 
face (or, more precisely, sight cell) which is selected based on its incidence, i.e. its ability to connect
edge crossing free to vertices incident to it. Hence, all subfaces identified for that face/sight cell will 
not induce any edge crossings for those incident vertices. All other columns will be non-zero.
"""


# The Data -------------------------------------------------------------------------------------------------------------

# Induced Edge Crossings of n_a = 7 subfaces (in face A) when connecting to the m = 4 target vertex neighbors
induced_edge_crossings_a = np.array(
    [
        [0, 1, 1, 1],
        [0, 1, 3, 4],
        [0, 3, 2, 2],
        [0, 4, 3, 1],
        [0, 1, 4, 1],
        [0, 2, 1, 4],
        [0, 1, 1, 3]
    ]
)

induced_edge_crossings_a_weights = np.random.uniform(low=0.0, high=5.0,
                                                     size=(induced_edge_crossings_a.shape[0],
                                                           induced_edge_crossings_a.shape[1]))

print(induced_edge_crossings_a_weights)

# Induced Edge Crossings of n_a = 4 subfaces (in face A) when connecting to the m = 4 target vertex neighbors
induced_edge_crossings_b = np.array(
    [
        [2, 0, 1, 1],
        [4, 0, 3, 2],
        [2, 0, 1, 1],
        [1, 0, 4, 1],
    ]
)

induced_edge_crossings_b_weights = np.random.uniform(low=0.0, high=5.0,
                                                     size=(induced_edge_crossings_b.shape[0],
                                                           induced_edge_crossings_b.shape[1]))
print(induced_edge_crossings_b_weights)
