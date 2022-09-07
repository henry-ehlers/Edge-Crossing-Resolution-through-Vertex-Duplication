import matplotlib.pyplot as plt
import networkx as nx
import itertools as it
import numpy as np
import timeit
import sys

from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
from src.diagnostics import *
from src.line_segments import *
from src.vertex_splitting import *
from src.sight_cells import *
from src.faces import *

# Specify vertices and edges
# todo: the example below causes floating point crashes as all their x and y points are identical
#coordinates = [(0, 0), (1, 2), (2, 0), (3, 2), (4, 0), (5, 3), (4, 1), (3, 3), (2, 1), (1, 3)]
coordinates = [(0, 2), (1, 0), (2, 1), (3, 0), (4, 2), (2, 4)]
target_vertices = [0, 2, 7]
vertices = range(0, len(coordinates))
edges = ((index, (index + 1) % len(vertices)) for index in range(0, len(vertices)))

more_coordinates = [(-2, 1.5), (-1, 1.5), (-1, 0.5)]
more_vertices = range(len(coordinates), len(coordinates) + len(more_coordinates))

more_edges = ((more_vertices[index], more_vertices[(index + 1) % len(more_vertices)])
              for index in range(0, len(more_vertices)))

# Create Graph
graph = nx.Graph()
for vertex in vertices:
    graph.add_node(vertex, target=1 if vertex in target_vertices else 0)
for vertex in more_vertices:
    print(vertex)
    graph.add_node(vertex, target=1 if vertex in target_vertices else 0)

for edge in edges:
    graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
for edge in more_edges:
    graph.add_edge(u_of_edge=edge[0], v_of_edge=edge[1], real=1)
positions = {vertices[index]: np.array(coordinates[index]) for index in range(0, len(coordinates))}
positions.update({more_vertices[index]: np.array(more_coordinates[index]) for index in range(0, len(more_vertices))})

# Create Output Directory
output_directory = "./drawings/tests/complex_outer_face"
Path(output_directory).mkdir(parents=True, exist_ok=True)

# Draw Initial Embedding
draw_graph(graph=graph, positions=positions)
save_drawn_graph(f"{output_directory}/sight_cell_line_segments_0.png")

# Planarize Graph
edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
plane_graph, plane_positions, virtual_edge_set = planarize_graph(
    graph=graph, positions=positions, edge_crossings=edge_crossings)

# # Locate faces and best two for target face
# TODO: find outer face
faces = find_all_faces(plane_graph, plane_positions)
face_edge_map = build_face_to_edge_map(plane_graph, faces)
face_incidences = find_face_vertex_incidence(faces, target_vertices)
ordered_face_edges = get_ordered_face_edges(faces, plane_graph)

# Outer Face
outer_faces = find_outer_face(ordered_face_edges, graph)
outer_face_identifier = frozenset(set.union(*[set(outer_face) for outer_face in outer_faces]))
outer_face_sorted_edges = [get_face_vertex_sequence(outer_face, plane_graph) for outer_face in outer_faces]
outer_face_sorted_vertices = [get_sorted_face_vertices(edge, is_sorted=True) for edge in outer_face_sorted_edges]

outer_cells, outer_edge_map = get_outer_face_sight_cells(selected_faces=outer_faces,
                                                         ordered_face_edges=ordered_face_edges,
                                                         graph=plane_graph,
                                                         positions=plane_positions,
                                                         bounds=((-6, -6), (-6, 6), (6, 6), (6, -6)))

# Draw and Save Planar rGraph
draw_graph(graph=plane_graph, positions=plane_positions)
save_drawn_graph(f"{output_directory}/sight_cell_line_segments_1.5.png")

outer_sight_cell_incidences = get_outer_face_sight_cell_incidences(sight_cells=outer_cells,
                                                                   target_vertices=target_vertices,
                                                                   face_edges=ordered_face_edges,
                                                                   face_edge_map=outer_edge_map,
                                                                   positions=plane_positions)

outer_sight_cell_edges = get_sight_cell_edges(outer_cells, plane_graph)

merge_cells_wrapper(face_sight_cells=outer_cells,
                    cells_edge_list=outer_sight_cell_edges,
                    cell_incidences=outer_sight_cell_incidences,
                    graph=plane_graph)

# Draw and Save Planar rGraph
draw_graph(graph=plane_graph, positions=plane_positions)
save_drawn_graph(f"{output_directory}/sight_cell_line_segments_1.75.png")