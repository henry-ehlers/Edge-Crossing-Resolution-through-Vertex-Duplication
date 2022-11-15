# Watts-Strogatz

## N = 10, E = 6, P = 0.5

- 5 (2 splits of nodes 8 and 9)
	- Adjacency of **Node 9**: 7, 3, 2, 5, 6, 4
	- Common Neighbors of Nodes **9** and **0**: 4, 6, 2
	- Are these Paths Possible
		- 1 > 3 > 9 > 4 [yes]
		- 2 > 3 > 4 > 8 [yay]
		- 8 > 7 > 0 > 9 [nay]
		- 9 > 5 > 8 > 0 [nay]

- 8 (2 splits: 3 and 7)
	- Adjacency of **Node 3**: 4, 1, 2, 0, 7, 8
	- Common Neighbors of Nodes **3** and **4**: 0, 1, 2
	- Are these Paths possible
		- 1 > 0 > 3 > 4 [yay]
		- 2 > 4 > 3 > 7 [yay]
		- 9 > 6 > 8 > 3 [nay]
		- 3 > 4 > 5 > 2 [nay]

- 39 (2 splits: 0 9)
	- Adjacency of Node **9**: 0, 6, 5, 7, 8, 3, 4
	- Common Neighbors of Nodes **0** and **9**: 5, 6, 3, 8, 7
	- Are these Paths possible:
		- 0 > 9 > 3 > 2 [yes]
		- 5 > 9 > 4 > 1 [yes]
		- 7 > 0 > 4 > 9 [no]
		- 1 > 8 > 9 > 5 [no]
		- 6 > 5 > 9 > 3 [yes]