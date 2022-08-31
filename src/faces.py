from src.edges import *
from src.line_segments import *
from src.edge_crossings import *
from src.vertex_splitting import calculate_face_centroid
from src.graph_drawing import *
from src.graph_drawing import get_graph_entity_data

import matplotlib.path as mpltPath
import scipy as sp
import networkx as nx
import itertools as it
import numpy as np
import math
import copy
import sys


def find_all_subfaces(graph, virtual_edge_set_map, target_face_to_vertex_map):

    # Prepare dictionary of sets within which to store all found faces per face
    subfaces = {target_face: set() for target_face in target_face_to_vertex_map.keys()}

    # Limit search of cycle basis to subfaces
    for target_face in target_face_to_vertex_map.keys():
        vertices_to_keep = set()

        # Mark vertices based on their inclusion in virtual edge sets pertaining to the current set
        for edge_set in virtual_edge_set_map:
            intersection = target_face_to_vertex_map[target_face] & edge_set
            if len(intersection) >= 2:
                [vertices_to_keep.add(vertex) for vertex in list(edge_set)]
            elif len(intersection) == 1:
                vertices_to_keep.add(intersection.pop())
            else:
                continue

        # Keep only vertices pertaining to current subgraph and search for cycle basis within
        subgraph = copy.deepcopy(graph).subgraph(list(vertices_to_keep))
        subfaces[target_face] = find_all_faces(subgraph)  # todo: also pass positions?

    # Return all subfaces for each target face
    return subfaces


def get_sight_cells_edge_sets(sight_cells, graph):
    sight_cell_edge_list = {face: {} for face in sight_cells.keys()}
    for face in sight_cells.keys():
        sight_cell_edge_list[face].update(get_sight_cell_edges(sight_cells[face], graph))
    return sight_cell_edge_list


def get_sight_cell_edges(sight_cells, graph):
    edge_set = {cell: set() for cell in sight_cells}
    for cell in sight_cells:
        cell_edges = get_face_vertex_sequence(cell, graph)
        [edge_set[cell].add(frozenset(edge)) for edge in cell_edges]
    return edge_set


def get_ordered_face_edges(faces, graph):
    ordered_face_edges = dict.fromkeys(faces)
    for face in faces:
        print(f"Face: {face}")
        ordered_face_edges[face] = get_face_vertex_sequence(face, graph)
    return ordered_face_edges


def shrink_cycle(cycle, other_cycles, sorted_edges, graph, positions):

    sorted_vertices = get_sorted_face_vertices(sorted_edges[cycle], is_sorted=True)
    cycle_coordinates = [positions[vertex] for vertex in sorted_vertices]

    # Get positions from polygon of ordered points of current cycle
    cycle_path = mpltPath.Path(cycle_coordinates[0:-1])
    for other_cycle in other_cycles:

        # Skip if the cycle == the other cycle
        if other_cycle == cycle:
            continue

        vertex_intersection = cycle.intersection(other_cycle)
        if len(vertex_intersection) == 0:
            continue

        remaining_vertices = other_cycle - vertex_intersection
        if len(remaining_vertices) == 0:
            continue

        remaining_coordinates = [positions[vertex] for vertex in remaining_vertices]
        in_side = cycle_path.contains_points(remaining_coordinates)
        if not all(in_side):
            if not all(not in_side[index] for index in range(0, len(remaining_vertices))):
                print("shit's fucked")
            continue

        #
        new_cycle, new_edge_list = None, None
        if len(vertex_intersection) == 1:
            # TODO: fix this section -> edge order is unclear
            new_cycle = cycle.union(remaining_vertices)
            new_edge_list = get_face_vertex_sequence(new_cycle, graph)

        elif len(vertex_intersection) >= 2:
            cycle_edge_list = set([frozenset(edge) for edge in sorted_edges[cycle]])
            other_edge_list = set([frozenset(edge) for edge in sorted_edges[other_cycle]])

            new_edge_set = cycle_edge_list.symmetric_difference(other_edge_list)
            new_edge_list = sort_face_edges(list([tuple(edge) for edge in new_edge_set]))
            new_cycle = frozenset().union(*new_edge_set)

        # Mention that something is broken
        if (new_cycle is None) or (new_edge_list is None):
            sys.exit("SHIT'S FUCKED")

        # Update the list of cycles
        other_cycles.remove(cycle)
        other_cycles.add(new_cycle)

        # Update the cycle sorted edge dictionary
        sorted_edges[new_cycle] = new_edge_list
        sorted_edges.pop(cycle)

        # Recurse and return
        return shrink_cycle(new_cycle, other_cycles, sorted_edges, graph, positions)


def find_singleton_cycles(cycles, graph, as_set=True):

    # Get Ordered Face edges for all identified cycles
    ordered_face_edges = get_ordered_face_edges(cycles, graph)

    # Find edges not Featured in Cycles
    non_cycle_edges = set()
    for edge in list(graph.edges):
        in_cycle = [set(edge) in [set(cycle_edge) for cycle_edge in ordered_face_edges[cycle]] for cycle in cycles]
        if any(in_cycle):
            continue
        non_cycle_edges.add(frozenset(edge))

    # Add found singleton cycles to identified ones
    if as_set:
        [cycles.add(singleton_cycle) for singleton_cycle in non_cycle_edges]
    else:
        [cycles.append(singleton_cycle) for singleton_cycle in non_cycle_edges]


def find_all_faces(graph, positions=None, as_set=True):

    # TODO: in all subsequent functions which use the face dictionary, check the length of the face (no singletons)

    # Identify the minimum cycle basis of the graph
    cycles = nx.minimum_cycle_basis(G=graph)
    cycles = set([frozenset(cycle) for cycle in cycles]) if as_set else [set(cycle) for cycle in cycles]

    # Check for and add Singleton edges as singleton cycles
    find_singleton_cycles(cycles, graph, as_set)

    # If no positions are given with which to check the validity of the cycles, return the cycles as faces
    if positions is None:
        return cycles

    # Get the ordered edge list per cycle
    print("a")
    cycle_ordered_edges = get_ordered_face_edges(cycles, graph)
    print("b")
    print(f"cycle ordered edges: {cycle_ordered_edges}")
    # For each Cycle recursively check that the face is indeed a face
    faces = copy.copy(cycles)
    for cycle in cycles:
        shrink_cycle(cycle=cycle,
                     other_cycles=faces,
                     sorted_edges=cycle_ordered_edges,
                     graph=graph,
                     positions=positions)

    # Return set of faces (frozen sets of vertices)
    return faces


def color_selected_faces(graph, face_set, face_edge_map):
    for face in face_set:
        for edge in face_edge_map[frozenset(face)]:
            graph.edges[edge]["target"] = 1
    return graph


def close_closed_face(convex_face_edge_list, graph):

    # Note down the face's edges
    face_edges = []

    # Iterate over all ORDERED vertices that make up the new convex face
    for index in range(0, len(convex_face_edge_list)):

        # Extract two adjacent vertices from the face
        vertex_a = convex_face_edge_list[index]
        vertex_b = convex_face_edge_list[(index + 1) % len(convex_face_edge_list)]

        # Check if an edge already exists between these vertices
        try:
            edge = graph.edges[vertex_a, vertex_b]
        except KeyError:
            edge = None

        # If no edge exists, add one and update the convex face's ordered edge list
        if edge is None:
            graph.add_edge(vertex_a, vertex_b, virtual=1, segment=0, target=0)

        # Append (new) edge to this face's list of edges
        face_edges.append((vertex_a, vertex_b))

    # Return update list of unsorted face edges
    return face_edges


def find_outer_face(ordered_face_edges, graph):

    # Initialize a dictionary which maps edges to the number of faces they are in
    faces_per_edge = dict.fromkeys(list(graph.edges()), 0)

    # Iterate over all faces and increment counts of edges within them
    for face in ordered_face_edges.keys():
        for edge in ordered_face_edges[face]:
            edge_a, edge_b = (edge[0], edge[1]), (edge[1], edge[0])
            if edge_a in faces_per_edge.keys():
                faces_per_edge[edge_a] += 1
            elif edge_b in faces_per_edge.keys():
                faces_per_edge[edge_b] += 1
            else:
                sys.exit("Shit's fucked")

    # Find edges which only map to a single face
    # min_face_count = min(faces_per_edge.values())
    faces = set([frozenset(edge) for edge in faces_per_edge.keys() if faces_per_edge[edge] == 1])
    print(f"singleton edges: {faces}")
    find_vertex_sets_from_edges(faces)

    # Return unique vertex sets from the found singleton edges
    return faces


def find_vertex_sets_from_edges(edge_sets):
    print(f"\nedge sets: {edge_sets}")
    for edge_set_a in edge_sets:
        print(f"Edge Set A: {edge_set_a}")
        for edge_set_b in edge_sets:
            print(f"Edge Set b: {edge_set_b}")
            if edge_set_a == edge_set_b:
                print("Continue")
                continue
            if not edge_set_a.intersection(edge_set_b):
                continue

            print("Merge")
            new_edge_set = edge_set_a.union(edge_set_b)
            print(f"new edge set: {new_edge_set}")
            edge_sets.remove(edge_set_a)
            edge_sets.remove(edge_set_b)
            edge_sets.add(new_edge_set)

            print(f"post merge sets: {edge_sets}")
            return find_vertex_sets_from_edges(edge_sets)



def build_face_to_edge_map(graph, faces):

    # TODO: the face set now consists of shortest cycles, so this basic look-up should be acceptable for now
    # TODO: however, for outer faces, there may be cycles within the the defined cycle that should not be included
    face_edge_map = dict()
    for face in faces:
        face_edge_map[face] = list()
        face_list = list(it.combinations(list(face), 2))
        for vertex_combination in face_list:
            vertex_a = vertex_combination[0]
            vertex_b = vertex_combination[1]
            found_edge = get_edge_connecting_vertices(graph, vertex_a, vertex_b)
            if found_edge: face_edge_map[face].append(found_edge)

    return face_edge_map


def get_edge_connecting_vertices(graph, vertex_a, vertex_b):
    for edge in graph.edges:
        if (vertex_a == edge[0] and vertex_b == edge[1]) or (vertex_b == edge[0] and vertex_a == edge[1]):
            return edge


def find_face_vertex_incidence(faces, target_vertices):

    # Initialize an empty dictionary (using sets as keys) to store vertex incidence sets
    face_incidences = dict()  # face_incidences = {face: dict() for face in faces}

    # Initialize set of vertex set as list, and target vertices as set
    faces = list([list(face) for face in faces])
    target_vertex_set = set(target_vertices)

    # Iterate over all faces
    for face_index_a in range(0, len(faces) - 1):

        # Extract current face, transform back into set
        face_a = frozenset(faces[face_index_a])
        face_incidences[face_a] = dict()

        # Determine how many target vertices are incident and left over
        incidence_a = face_a & target_vertex_set
        remaining_targets = target_vertex_set - incidence_a

        # Iterate over remaining faces
        for face_index_b in range(face_index_a + 1, len(faces)):

            # Extract current face, transform back into set
            face_b = frozenset(faces[face_index_b])

            # Determine how many remaining target faces are incident
            incidence_b = face_b & remaining_targets
            incident_vertices = incidence_a.union(incidence_b)

            # Store results only if non-empty
            if incident_vertices:
                face_incidences[face_a][face_b] = incident_vertices

    # Return Dictionary of face vertex incidence
    return face_incidences


def get_maximally_incident_faces(face_incidences):
    max_incidence = float("-inf")
    selected_faces = []
    for face_a in face_incidences.keys():
        for face_b in face_incidences[face_a].keys():
            vertex_incidence = len(face_incidences[face_a][face_b])
            if vertex_incidence > max_incidence:
                max_incidence = vertex_incidence
                selected_faces = [[face_a, face_b]]
            elif vertex_incidence == max_incidence:
                selected_faces.append([face_a, face_b])

    return max_incidence, selected_faces


def cross_product(vector_a, vector_b):
    return (vector_a[0] * vector_b[1]) - (vector_b[0] * vector_a[1])


def vector_angle(vector_1, vector_2):

    # Calculate Unit Vectors of Input Vectors
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)

    # Calculate Dot Product and Signed Angle in Radians
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(dot_product)

    # Return Angle in Degrees
    return np.degrees(angle)


# Calculate Signed Area of Polygon
def signed_area(ordered_points_list):
    x, y = [point[0] for point in ordered_points_list], [point[1] for point in ordered_points_list]
    return sum(x[i] * (y[i + 1] - y[i - 1]) for i in range(-1, len(x) - 1)) / 2.0


def calculate_face_signed_area(ordered_face_vertices, positions):
    vertex_points = [positions[vertex] for vertex in ordered_face_vertices]
    area = signed_area(vertex_points)
    return area


def calculate_outer_angle(point_a, point_b, point_c):

    # Assumed that point b connects to both a and c
    vector_1, vector_2 = point_a - point_b, point_c - point_b

    # Calculate Signed Angle between two Vectors
    signed_angle = vector_angle(vector_1, vector_2)

    # Return Outer angle
    return signed_angle if cross_product(vector_1, vector_2) < 0.0 else 360 - signed_angle


def calculate_face_outer_angles(counter_clockwise_face_vertices, positions):

    # Initialize emtpy dictionary to store inner angles
    outer_angles = {vertex: None for vertex in counter_clockwise_face_vertices}

    # Iterate over all vertices, and calculate angle for each
    for vertex_index in range(0, len(counter_clockwise_face_vertices)):

        # Get three vertices that form an angle (the ordered hereof is crucial)
        vertex_a = counter_clockwise_face_vertices[vertex_index]
        vertex_b = counter_clockwise_face_vertices[vertex_index - 1]
        vertex_c = counter_clockwise_face_vertices[vertex_index - 2]

        # Calculate Inner angle and store as with center vertex as key
        point_a, point_b, point_c = positions[vertex_a], positions[vertex_b], positions[vertex_c]
        outer_angles[vertex_b] = calculate_outer_angle(point_a, point_b, point_c)

    # Return all inner angles
    return outer_angles


def calculate_inner_angle(point_a, point_b, point_c):

    # Assumed that point b connects to both a and c
    vector_1, vector_2 = point_a - point_b, point_c - point_b

    # Calculate Signed Angle between two Vectors
    signed_angle = vector_angle(vector_1, vector_2)

    # Return inner angle
    return signed_angle if cross_product(vector_1, vector_2) > 0.0 else 360 - signed_angle


def calculate_face_inner_angles(counter_clockwise_face_vertices, positions):

    # Initialize emtpy dictionary to store inner angles
    inner_angles = {vertex: None for vertex in counter_clockwise_face_vertices}

    # Iterate over all vertices, and calculate angle for each
    for vertex_index in range(0, len(counter_clockwise_face_vertices)):

        # Get three vertices that form an angle (the ordered hereof is crucial)
        vertex_a = counter_clockwise_face_vertices[vertex_index]
        vertex_b = counter_clockwise_face_vertices[vertex_index - 1]
        vertex_c = counter_clockwise_face_vertices[vertex_index - 2]

        # Calculate Inner angle and store as with center vertex as key
        point_a, point_b, point_c = positions[vertex_a], positions[vertex_b], positions[vertex_c]
        inner_angles[vertex_b] = calculate_inner_angle(point_a, point_b, point_c)

    # Return all inner angles
    return inner_angles


def add_boundary_to_graph(bounds, graph, positions, offset=0.2):

    # Define the labels of vertices and their edges
    bound_vertices = list(range(max(graph.nodes) + 1, max(graph.nodes) + 1 + len(bounds)))
    bound_edges = [(bound_vertices[ind], bound_vertices[(ind + 1) % len(bounds)]) for ind in range(0, len(bounds))]

    # Update Graph, Edges, and Vertex Positions
    for index in range(0, len(bounds)):
        new_position = np.array(bounds[index]) - np.sign(np.array(bounds[index])) * np.array([offset, offset])
        positions[bound_vertices[index]] = new_position
    graph.add_nodes_from(bound_vertices)
    graph.add_edges_from(bound_edges)
    # Return the added vertex labels and edges
    return bound_vertices, bound_edges


def get_face_sight_cells(selected_faces, ordered_face_edges, graph, positions,
                         bounds=((-1, -1), (-1, 1), (1, 1), (1, -1)),
                         outer=False):

    bound_vertices, bound_edges = add_boundary_to_graph(bounds, graph, positions) if outer else ([], [])


    #
    face_edge_map = {face: {} for face in selected_faces}
    sight_cells = {face: None for face in selected_faces}

    # Iterate over all faces
    for face in selected_faces:

        # Add additional candidate edges if we are dealing with the outer face
        other_faces = copy.copy(selected_faces)
        other_faces.remove(face)
        print(f"Others = {other_faces} = {selected_faces} - {face}")
        print(f"All edges: {ordered_face_edges}")
        other_face_edges = unlist([ordered_face_edges.get(other_face) for other_face in other_faces]) if outer else []
        print(f"For Face {face} we additionally consider edge {other_face_edges}")

        # Get Vertices and ensure they are listed in counter-clockwise order
        face_edges = ordered_face_edges[face]
        face_vertices = get_sorted_face_vertices(face_edges, is_sorted=True)
        if calculate_face_signed_area(face_vertices, positions) < 0:
            face_vertices = list(reversed(face_vertices))

        # Calculate Inner Angles to check convexity
        face_angles = calculate_face_inner_angles(face_vertices, positions) if not outer \
            else calculate_face_outer_angles(face_vertices, positions)
        print(f"Face Angles: {face_angles}")
        # If the face is convex, the face is the sight-cell
        if is_convex(face_angles) and not outer:
            face_edge_map[face].update({edge: [edge] for edge in face_edges})
            sight_cells[face] = face_edges

        # Otherwise, split into sight cells and return set of frozen sets of vertices per sight-cell
        else:
            candidate_edges = face_edges + bound_edges + other_face_edges
            print(f"Candidate Edges: {candidate_edges}")
            cells, virtual_edge_map = project_face_sight_lines(edges=candidate_edges,
                                                               vertices=face_vertices,
                                                               inner_angles=face_angles,
                                                               graph=graph,
                                                               positions=positions,
                                                               bounds=bounds,
                                                               outer=outer)
            face_edge_map[face].update(virtual_edge_map)
            sight_cells[face] = cells

            # Also add those edges that were not bisected
            face_edge_map[face].update({edge: [edge] for edge in face_edges if edge not in virtual_edge_map.keys()})

    # Return
    return sight_cells, face_edge_map


def get_outer_face_sight_cells(selected_faces, ordered_face_edges, graph, positions,
                         bounds=((-1, -1), (-1, 1), (1, 1), (1, -1))):

    bound_vertices, bound_edges = add_boundary_to_graph(bounds, graph, positions)
    all_face_edges = unlist([ordered_face_edges.get(face) for face in selected_faces])

    #
    face_edge_map = {edge: [edge] for edge in all_face_edges + bound_edges}
    sight_cells = {}
    # TODO: CANDIDATE EDGES ARE ALWAYS THE SAME

    # Iterate over all faces
    for face in selected_faces:

        # Add additional candidate edges if we are dealing with the outer face
        other_faces = copy.copy(selected_faces)
        other_faces.remove(face)


        # Get Vertices and ensure they are listed in counter-clockwise order
        face_edges = unlist([face_edge_map.get(edge) for edge in ordered_face_edges[face]])
        print(f"face edges: {face_edges}")
        face_vertices = get_sorted_face_vertices(face_edges, is_sorted=True)
        if calculate_face_signed_area(face_vertices, positions) < 0:
            face_vertices = list(reversed(face_vertices))

        # Calculate Inner Angles to check convexity
        face_angles = calculate_face_outer_angles(face_vertices, positions)
        print(f"face angles: {face_angles}")
        print(f"candidate edges: {all_face_edges + bound_edges}")
        print(f"map: {face_edge_map}")
        candidate_edges = unlist([face_edge_map.get(edge) for edge in all_face_edges + bound_edges])
        print(f"Candidate Edges: {candidate_edges}")
        cells, virtual_edge_map = project_face_sight_lines(edges=candidate_edges,
                                                           vertices=face_vertices,
                                                           inner_angles=face_angles,
                                                           graph=graph,
                                                           positions=positions,
                                                           bounds=bounds,
                                                           outer=True)
        print(f"virtual edge map: {virtual_edge_map}")
        face_edge_map.update(virtual_edge_map)

        # Also add those edges that were not bisected
        face_edge_map.update({edge: [edge] for edge in face_edges + bound_edges
                              if edge not in virtual_edge_map.keys()})

    print(f"face map: {face_edge_map}")
    # Find Sight Cells
    # define subset of graph to be investigated
    outer_face_vertices = list(set(unlist(unlist(list(face_edge_map.values())))))
    print(f"outer face vertices: {outer_face_vertices}")
    face_graph = nx.Graph(graph.subgraph(nodes=outer_face_vertices))
    print(f"\n shits fucked")
    print(f"vertices: {face_graph.nodes}")
    print(f"edges: {face_graph.edges}")
    sight_cells = find_all_faces(face_graph, positions=positions)

    # remove map-corrected vertex sets corresponding to original outer faces from cells
    print(f"\nsight cells: {sight_cells}")
    sys.exit()
    # Return
    return sight_cells, face_edge_map


def is_convex(inner_angles):
    return all(angle <= 180.0 for angle in inner_angles.values())


def are_vertices_adjacent(vertex_a, vertex_b, graph):
    adjacent = True
    try:
        graph.edges[vertex_a, vertex_b]
    except KeyError:
        adjacent = False
    return adjacent


def project_face_sight_lines(edges, vertices, inner_angles, graph, positions, bounds, outer):

    # Keep track of the added vertices, and in which edges they were added
    added_vertices, edge_to_virtual_vertices = [], {}
    # Consider only those vertices whose angle is greater than 180 degrees
    bend_vertices = [key for key in inner_angles.keys() if inner_angles[key] > 180]

    for joint_vertex in bend_vertices:
        print(f"\njoint vertex: {joint_vertex}")
        for connecting_vertex in vertices:
            print(f"connecting vertex: {connecting_vertex}")

            # Skip any vertex pair that is a) consists of the same vertex, or b) has already been investigated
            if connecting_vertex == joint_vertex:
                continue

            # Check whether bend and other vertex can 'see' each other
            is_visible = True if are_vertices_adjacent(joint_vertex, connecting_vertex, graph) \
                else is_vertex_visible(joint_vertex=joint_vertex,
                                       connecting_vertex=connecting_vertex,
                                       inner_angles=inner_angles,
                                       graph=graph,
                                       vertices=vertices,
                                       edges=edges,
                                       positions=positions,
                                       outer=outer)

            # If they cannot see each other, skip to the next pair
            if not is_visible:
                continue
            print("is visible")
            # Extend the sight-line, producing a
            # TODO: extend "outer" keyword all the way down
            bisected_edge, new_vertex = extend_sight_line(joint_vertex=joint_vertex,
                                                          connecting_vertex=connecting_vertex,
                                                          inner_angles=inner_angles,
                                                          vertices=vertices,
                                                          edges=edges,
                                                          graph=graph,
                                                          positions=positions,
                                                          bounds=bounds,
                                                          outer=outer)
            print(f"new vertex: {new_vertex}")
            # Keep track of what has been added
            added_vertices.append(new_vertex)

            # Vertices can see one-another, but not produce a legal extension.
            if bisected_edge is None:
                continue

            # Add
            print(f"Bisected Edge: {bisected_edge}")
            if bisected_edge in edge_to_virtual_vertices:
                edge_to_virtual_vertices[bisected_edge].add(new_vertex)
            else:
                edge_to_virtual_vertices[bisected_edge] = {new_vertex}

    print(f"edge to virtual: {edge_to_virtual_vertices}")

    # Remove edges which have been intersected, and replace them with ordered virtual edges
    virtual_edge_map = add_virtual_edges(graph, positions, edge_to_virtual_vertices)
    print(f"virtual edge map: {virtual_edge_map}")
    remove_edges(graph, edge_to_virtual_vertices.keys())

    # Locate Edge Crossings and Faces in Subgraph
    candidate_edge_vertices = list(set(unlist(edges)))
    print(f"added vertices: {added_vertices}")
    print(f"candidate edge vertices: {candidate_edge_vertices}")
    face_vertices = vertices + added_vertices + candidate_edge_vertices
    face_positions = {key: positions.get(key) for key in face_vertices}
    face_graph = nx.Graph(graph.subgraph(nodes=face_vertices))

    # Find remaining edge crossings (between placed line-segments) and replace them with virtual vertices
    face_edge_crossings, vertex_crossings = locate_edge_crossings(face_graph, face_positions)
    print(f"edge crossings: {face_edge_crossings}")
    if len(face_edge_crossings) > 1:
        face_graph, face_positions, virtual_edges = planarize_graph(face_graph, face_positions, face_edge_crossings)
        non_empty_virtual_edges = {k: v for k, v in virtual_edges.items() if v}
        print(f"removed: {non_empty_virtual_edges}")
        virtual_edge_map.update(non_empty_virtual_edges)
        print(f"virtual edges: {virtual_edges}")
        graph.update(face_graph)
        positions.update(face_positions)
        [graph.remove_edge(u=edge[0], v=edge[1]) for edge in virtual_edges.keys() if virtual_edges[edge]]

    draw_graph(graph=graph, positions=positions)
    save_drawn_graph("sight_cell_line_segments_1.25.png")

    # Define Sight Cells, i.e. faces
    sight_cells = find_all_faces(face_graph, positions=positions) if not outer else None

    # Return Sight Cells
    return sight_cells, virtual_edge_map


def unlist(nested_list):
    return list(it.chain.from_iterable(nested_list))


def get_face_sight_cell_incidences(sight_cells, face_incidences, target_vertices, face_edges, face_edge_map, positions):

    # Initialize an empty map of sight cell to incidence
    sight_cell_incidences = {sight_cell: set() for sight_cell in list(it.chain.from_iterable(sight_cells.values()))}

    # Iterate over all faces in the graph
    for face in sight_cells.keys():

        # If the original face is convex, the sight cell is equivalent
        if len(sight_cells[face]) == 1:
            sight_cell_incidences[face].add(face_incidences[face])
            continue

        # Extract all edges in the face, i.e. the virtual edges formed by virtual edge bisection
        face_edge_list = unlist([face_edge_map[face][edge] for edge in face_edges[face]])

        # Iterate over very sight cell in the current face and check visibility to the target vertices
        for sight_cell in sight_cells[face]:
            sight_cell_incidences[sight_cell].update(get_sight_cell_incidence(sight_cell_vertices=sight_cell,
                                                                              target_vertices=target_vertices,
                                                                              real_face_edges=face_edge_list,
                                                                              positions=positions))

    # Return a
    return sight_cell_incidences


def get_sight_cell_incidence(sight_cell_vertices, target_vertices, real_face_edges, positions):

    # Check what targets are already in the current sight cell
    targets_in_cell = sight_cell_vertices.intersection(target_vertices)
    remaining_targets = set(target_vertices) - targets_in_cell

    # Define the sight cell's positions and centroid, as well as initialize incidence set
    sight_cell_positions = [positions[vertex] for vertex in list(sight_cell_vertices)]
    sight_cell_centroid = calculate_face_centroid(sight_cell_positions)
    sight_cell_incidence = set(targets_in_cell)

    # Iterate over all targets that are not already incident to the cell
    for target in remaining_targets:
        sight_line = [sight_cell_centroid, positions[target]]

        # Iterate over all edges and check whether they intersect the line between centroid and target
        is_target_visible = True
        for edge in real_face_edges:
            edge_line = [positions[vertex] for vertex in edge]
            intersection = line_intersection(sight_line[0], sight_line[1], edge_line[0], edge_line[1])

            # If the two lines intersect, they are not visible
            if intersection is not None:
                is_target_visible = False
                break

        # If no edges intersected the sight line, append it as a visible target
        if is_target_visible:
            sight_cell_incidence.add(target)

    # Return the found set of visible vertices
    return sight_cell_incidence


def merge_face_sight_cells(cells, cells_edge_list, cell_incidences, removed_vertices, graph):
    # todo: outer face also contains all faces INSIDE of the graph
    # iterate over all pairs of sight cells
    # for each pair check if incidence is the same, and they have an edge in common
    # if so, merge the two in the "sight_cells" object and delete the common edge
    # also delete one of the vertices that formed the edge in question
    # recurse back into the same function if at least one thing was merged, else return sight cells

    # Iterate over all cell pairs
    for cell_index_a in range(0, len(cells)):
        for cell_index_b in range(cell_index_a + 1, len(cells)):
            cell_a, cell_b = cells[cell_index_a], cells[cell_index_b]

            # Attempt to merge the two cells, return a boolean for success and a (possibly empty) vertex ID
            merge_successful, removed = try_merge_two_sight_cells(cell_a=cell_a,
                                                                           cell_b=cell_b,
                                                                           cells=cells,
                                                                           cells_edge_list=cells_edge_list,
                                                                           cell_incidences=cell_incidences,
                                                                           graph=graph)

            # If the merge was successful, recurse
            if merge_successful:

                # If a vertex was removed, keep track if its removal
                if removed:
                    [removed_vertices.append(removed_vertex) for removed_vertex in removed]

                # Recurse and repeat merging
                removed_vertices = merge_face_sight_cells(cells=cells,
                                                          cells_edge_list=cells_edge_list,
                                                          cell_incidences=cell_incidences,
                                                          removed_vertices=removed_vertices,
                                                          graph=graph)

                # Exit recursion (as the set of cells has changed) and return the removed vertices
                return removed_vertices

    # Final loop, nothing left to merge, return the removed vertices
    return removed_vertices


def try_merge_two_sight_cells(cell_a, cell_b, cells, cells_edge_list, cell_incidences, graph):

    # Determine along which the two cells are to be merged and where the cells are located in the list
    merge_edges = cells_edge_list[cell_a].intersection(cells_edge_list[cell_b])
    incidence_a, incidence_b = cell_incidences[cell_a], cell_incidences[cell_b]
    non_overlapping_incidences = incidence_a ^ incidence_b

    if non_overlapping_incidences or len(merge_edges) == 0:
        return False, None

    # Determine the new cell's vertex set and edge list
    new_edge_set = cells_edge_list[cell_a].union(cells_edge_list[cell_b] - merge_edges)
    new_cell = cell_a.union(cell_b)

    # Update the Cell List
    cells.append(new_cell)
    [cells.remove(cell_key) for cell_key in [cell_a, cell_b]]

    # Update each Cell's edge list
    cells_edge_list[new_cell] = new_edge_set
    [cells_edge_list.pop(cell_key, None) for cell_key in [cell_a, cell_b]]

    # Update the Cells Incidences
    cell_incidences[new_cell] = copy.deepcopy(cell_incidences[cell_a])
    [cell_incidences.pop(cell_key, None) for cell_key in [cell_a, cell_b]]

    # Update the graph
    common_vertices = []
    for merge_edge in merge_edges:
        update_merge_sight_cell_graph(merge_edge, graph)

    # The merge has successfully occurred
    return True, common_vertices


def update_merge_sight_cell_graph(merge_edge, graph):

    merge_edge = list(merge_edge)
    graph.remove_edge(u=merge_edge[0], v=merge_edge[1])

    # Identify a vertex connected to deleted edge which is couched between virtual edges. Remove vertex and replace edge
    vertex_edges, common_vertex = unlist([list(graph.edges(v)) for v in merge_edge if len(graph.edges(v)) == 2]), None

    # todo: maybe clean up the graph of non-connected and virtual edges?
    # if vertex_edges:
    #     # Define a new edge which skips the now singleton vertex
    #     new_edge = [vertex for vertex in unlist(vertex_edges) if vertex not in merge_edge]
    #     real_edge = all([get_graph_entity_data(graph.nodes, vertex, "virtual", 0) for vertex in new_edge])
    #     graph.add_edge(u_of_edge=new_edge[0], v_of_edge=new_edge[1], virtual=0 if real_edge else 1)
    #     print(f"Removed Edges {vertex_edges} with {new_edge}")
    #     # Remove the singleton virtual vertex from the graph and cells
    #     common_vertex = set(vertex_edges[0]).intersection(set(vertex_edges[1])).pop()
    #     graph.remove_node(common_vertex)
    #     print(f"REMOVING VERTEX {common_vertex}")
    #     return common_vertex

    return None


def merge_all_face_cells(face_sight_cells, face_cell_edge_map, cell_incidences, graph):

    # Iterate over all faces, and attempt to merge their faces
    for face in face_sight_cells.keys():

        # Skip convex faces
        if len(face_sight_cells[face]) == 1:
            continue

        # Try Merging Cells in non-convex face
        removed_vertices = merge_face_sight_cells(cells=list(face_sight_cells[face]),
                                                  cells_edge_list=face_cell_edge_map[face],
                                                  cell_incidences=cell_incidences,
                                                  removed_vertices=[],
                                                  graph=graph)

        # Update the face's cells, their incidents, and edges based on deleted vertices
        update_merged_sight_cell_data(face_cells=face_sight_cells[face],
                                      face_cell_incidences=cell_incidences,
                                      face_cell_edges=face_cell_edge_map[face],
                                      deleted_vertices=frozenset(removed_vertices))


def update_merged_sight_cell_data(face_cells, face_cell_incidences, face_cell_edges, deleted_vertices):
    """

    :param face_cells:
    :param face_cell_incidences:
    :param face_cell_edges:
    :param deleted_vertices: a FROZENSET of vertex ID's
    :return: nothing; the provided three datastructure are modified in place and returned by reference
    """

    # Skip if no vertices were deleted
    if not deleted_vertices:
        return

    # Clear Face's Old set of sight cells
    face_cells.clear()

    # Update Face-Level Information
    old_cells = list(face_cell_incidences.keys())
    for old_cell in old_cells:

        # Check whether the old cell contains any deleted vertices
        if not deleted_vertices.intersection(old_cell):
            continue

        # Create new cell without any of the previously deleted vertices
        new_cell = set(old_cell)
        [new_cell.discard(vertex) for vertex in deleted_vertices]
        new_cell = frozenset(new_cell)

        # Update The List of Cells
        face_cells.add(new_cell)

        # Update the Face Cell Incidences
        face_cell_incidences[new_cell] = face_cell_incidences.pop(old_cell)

        # Update the Dictionary of Cell Edges
        face_cell_edges[new_cell] = face_cell_edges.pop(old_cell)
        edges = copy.copy(face_cell_edges[new_cell])
        [face_cell_edges[new_cell].remove(edge) for edge in edges if deleted_vertices.intersection(edge)]


def find_minimal_sight_cell_set(cell_incidences):
    """

    :param cell_incidences: SUBSET OF ALL CELL INCIDENCES FOR A PARTICULAR FACE
    :return:
    """
    # TODO: implement Anais' thing here instead. we don;t need one task per one worker, but multiple per worker
    # todo: get target vertex list from list of unique vertices in incidences
    # Initiate cost matrix of ones
    # cost_matrix = np.ones(shape=(len(face_cells_incidence), len(target_vertices)), dtype=int)
    #
    # # Store sight cells in list to avoid ordering problems
    # row_names = list(face_cells_incidence.keys())
    # # Iterate over all sight cells and extract their incidences to build cost matrix
    # for sight_cell in row_names:
    #     row_index = row_names.index(sight_cell)
    #     for visible_vertex in face_cells_incidence[sight_cell]:
    #         col_index = target_vertices.index(visible_vertex)
    #         cost_matrix[row_index, col_index] -= 1

    # Find minimal assignment cost
    incidence_number = {cell: len(cell_incidences[cell]) for cell in cell_incidences.keys()}
    incidence_number = dict(sorted(incidence_number.items(), key=lambda item: item[1], reverse=True))
    selected_cells = list(incidence_number)[0:2]  # todo: selection could be arbitrarily long
    return {cell: cell_incidences[cell] for cell in selected_cells}


def select_sight_cells(sight_cells, sight_cell_incidence):
    minimal_cell_set = {face: {} for face in sight_cells.keys()}
    for face in sight_cells.keys():
        face_cell_incidences = {cell: sight_cell_incidence[cell] for cell in sight_cells[face]}
        minimal_cell_set[face].update(find_minimal_sight_cell_set(face_cell_incidences))
    return minimal_cell_set


def match_cell_and_face_incidence(face_incidences, selected_sight_cell_incidences):

    # Indicator whether faces must be reranked (i.e if incidences did not match)
    rerank_faces = False

    # Iterate over all (both selected faces)
    for face in face_incidences.keys():

        # Extract the set of target vertices incident to the current face and its selected cells
        face_incidence, cell_incidences = face_incidences[face], selected_sight_cell_incidences[face]

        # If the number of cells is greater than 1, than no single face contains the target incidence set
        if len(cell_incidences) > 1:

            # Replace the original face with its sight cells
            face_incidences.pop(face)
            face_incidences.update(cell_incidences)

            # Move on to next face and rerank all faces.
            rerank_faces = True
            break

    # Return indicator
    return rerank_faces


def extend_sight_line(joint_vertex, connecting_vertex, inner_angles, vertices, edges, graph, positions, bounds, outer):

    # Calculate intersections of extended line with boundaries in both directions
    bound_intersections = extend_line(positions[joint_vertex], positions[connecting_vertex], bounds)
    closest_intersection_to_joint = bound_intersections[0]
    print(f"bound intersection: {bound_intersections}")

    # If vertices are adjacent, they can see one-another; otherwise we must check explicitly
    already_connected = are_vertices_adjacent(joint_vertex, connecting_vertex, graph)
    is_visible = True \
        if already_connected \
        else check_vertex_visibility_by_angle(joint_vertex=joint_vertex,
                                              inner_angles=inner_angles,
                                              edges=edges,
                                              vertices=vertices,
                                              positions=positions,
                                              connecting_position=closest_intersection_to_joint,
                                              outer=outer)

    # If the hypothetical and observed angle are incompatible, then continue
    if not is_visible:
        return None, None
    print("extension visible")\

    # Find the Closest Intersection of the extended line with edges not incident to joint or connecting vertex
    extended_line = (positions[joint_vertex], closest_intersection_to_joint)
    print(f"extended_line line: {extended_line}")
    candidate_edges = [edge for edge in edges if not set(edge).intersection((joint_vertex, connecting_vertex))]
    print(f"candidate edges: {candidate_edges}")
    print(f"positions: {positions}")
    closest_edge, crossing_point = find_closest_edge_intersection(extended_line, candidate_edges, positions)

    # Add Virtual Vertex at Point of Intersection and a virtual edge between it and the origin
    origin_vertex, new_vertex_index = joint_vertex, max(graph.nodes) + 1
    graph.add_node(node_for_adding=new_vertex_index, virtual=1)
    graph.add_edge(u_of_edge=origin_vertex, v_of_edge=new_vertex_index, segment=1)
    positions[new_vertex_index] = np.array(crossing_point)
    if not already_connected:
        graph.add_edge(u_of_edge=joint_vertex, v_of_edge=connecting_vertex, segment=1)

    # Return a list of added vertices and a map of edges to newly placed virtual vertices
    return closest_edge, new_vertex_index


def is_vertex_visible(joint_vertex, connecting_vertex, inner_angles, graph, vertices, edges, positions, outer):

    # If vertices are neighbors, they can see one-another
    if are_vertices_adjacent(joint_vertex, connecting_vertex, graph):
        return True

    # Check Angle around the Joint Vertex allows for visibility to the connecting vertex
    angle_visibility = check_vertex_visibility_by_angle(joint_vertex=joint_vertex,
                                                        inner_angles=inner_angles,
                                                        edges=edges,
                                                        vertices=vertices,
                                                        positions=positions,
                                                        connecting_vertex=connecting_vertex,
                                                        outer=outer)

    # If the angle does not allow for visibility, return False
    if not angle_visibility:
        return angle_visibility

    # If the angle (hypothetically) allows for visibility, now check all possible edges that could be in the way
    possible_crossing_edges = [edge for edge in edges if (joint_vertex not in edge) and (connecting_vertex not in edge)]
    crossing_visibility = check_vertex_visibility_by_crossing(vertex_a=joint_vertex,
                                                              vertex_b=connecting_vertex,
                                                              candidate_edges=possible_crossing_edges,
                                                              positions=positions)

    # If the line segment between joint and connecting vertex crosses and edge, return False
    if not crossing_visibility:
        return crossing_visibility

    # All Conditions met for Visibility, so return True
    return True


def check_vertex_visibility_by_angle(joint_vertex, inner_angles, edges, vertices, positions, outer,
                                     connecting_vertex=None, connecting_position=None):

    # Ensure that either a vertex or a position has been provided
    assert (connecting_vertex is not None) or (connecting_position is not None), \
        "Specify a Connecting Point when checking its visibility by angle"

    # Get points for new angle calculation
    joint_index = [index for index in range(0, len(vertices)) if vertices[index] == joint_vertex][0]
    ref_vertex_a, ref_vertex_b = vertices[(joint_index + 1) % len(vertices)], vertices[joint_index - 1]
    debug_angle = calculate_outer_angle(positions[ref_vertex_a], positions[joint_vertex], positions[ref_vertex_b]) \
        if outer else calculate_inner_angle(positions[ref_vertex_a], positions[joint_vertex], positions[ref_vertex_b])

    # Get the Angle of the Joint against which we are comparing the new incoming angle:
    observed_angle = inner_angles[joint_vertex]
    assert observed_angle == debug_angle, \
        "Angle Calculation is incorrect"

    # Calculate Hypothetical Angle
    connecting_position = connecting_position if (connecting_position is not None) else positions[connecting_vertex]
    hypothetical_angle = calculate_outer_angle(positions[ref_vertex_a], positions[joint_vertex], connecting_position) \
        if outer else calculate_inner_angle(positions[ref_vertex_a], positions[joint_vertex], connecting_position)
    # If the angle between Vertex A and B is larger than between Vertex A and its neighbors, Vertex B is not visible
    return hypothetical_angle < observed_angle


def check_vertex_visibility_by_crossing(vertex_a, vertex_b, candidate_edges, positions):

    # Get Points of Vertices A and B
    point_a, point_b = positions[vertex_a], positions[vertex_b]

    # For each candidate edge, check crossing between it and (Vertex A, Vertex B)
    for edge in candidate_edges:

        # If Lines intersect, Vertices A and B cannot see one-another
        point_c, point_d = positions[edge[0]], positions[edge[1]]
        if line_intersection(point_a, point_b, point_c, point_d) is not None:
            return False

    # If no real edges intersected
    return True
