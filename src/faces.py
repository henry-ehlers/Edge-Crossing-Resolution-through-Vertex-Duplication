from src.edges import *
from src.graph_drawing import *

import matplotlib.path as mpltPath
import networkx as nx
import itertools as it
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import re
import copy
import sys


def get_face_sub_face_edge_sets(face_sub_cells, graph):
    sight_cell_edge_list = {face: {} for face in face_sub_cells.keys()}
    for face in face_sub_cells.keys():
        print(f"\nface: {face}")
        sight_cell_edge_list[face].update(get_sub_face_edges(face_sub_cells[face], graph))
    return sight_cell_edge_list


def get_sub_face_edges(sub_faces, graph):
    edge_set = {sub_face: set() for sub_face in sub_faces}
    for sub_face in sub_faces:
        cell_edges = get_face_vertex_sequence(sub_face, graph)
        [edge_set[sub_face].add(frozenset(edge)) for edge in cell_edges]
    return edge_set


def sub_face_merge_wrapper(face_sub_faces, face_sub_face_crossings, face_sub_face_edge_set,
                           graph, positions):

    #face_sub_face_edge_map

    for face in face_sub_faces.keys():
        # Try Merging Cells in non-convex face
        sub_faces = list(face_sub_faces[face])  # Reminder: cast is done to enable more efficient indexed looping
        merge_sub_faces(sub_faces=sub_faces,
                        sub_face_edge_set=face_sub_face_edge_set[face],
                        sub_face_crossings=face_sub_face_crossings[face],
                        graph=graph)
        face_sub_faces[face] = set(sub_faces)

        # Draw Merged Embedding
        draw_graph(graph=graph, positions=positions)
        save_drawn_graph(f"./merged_sub_faces_{face}.png")

    # Return updated sight cells, incidences, and edge map
    return face_sub_faces


def merge_sub_faces(sub_faces, sub_face_edge_set, sub_face_crossings, graph):

    #
    for index_a in range(0, len(sub_faces) - 1):
        for index_b in range(index_a + 1, len(sub_faces)):
            sub_face_a, sub_face_b = sub_faces[index_a], sub_faces[index_b]

            # Attempt to merge the two cells, return a boolean for success and a (possibly empty) vertex ID
            merge_successful = try_merge_two_sub_faces(sub_face_a=sub_face_a,
                                                       sub_face_b=sub_face_b,
                                                       sub_faces=sub_faces,
                                                       sub_face_edge_set=sub_face_edge_set,
                                                       sub_face_crossings=sub_face_crossings,
                                                       graph=graph)

            # If the merge was successful, recurse
            if merge_successful:

                # Recurse and repeat merging
                merge_sub_faces(sub_faces=sub_faces,
                                sub_face_edge_set=sub_face_edge_set,
                                sub_face_crossings=sub_face_crossings,
                                graph=graph)

                # Exit recursion (as the set of cells has changed) and return the removed vertices
                return

    # Final loop, nothing left to merge, return the removed vertices
    return


def try_merge_two_sub_faces(sub_face_a, sub_face_b, sub_faces, sub_face_edge_set, sub_face_crossings, graph):
    """

    :param sub_face_a: frozenset(sub_face)
    :param sub_face_b:
    :param sub_faces: list of frozenset(sub_faces)
    :param sub_face_edge_set: dictionary which maps sub_faces to its edges. The dictionary is typed and structured as
    follows: {frozenset(sub_face: (frozenset(edges)}
    :param sub_face_crossings:
    :param graph:

    :return:
    """

    # Check whether the two sub_faces have at least one edge in common
    merge_edges = sub_face_edge_set[sub_face_a].intersection(sub_face_edge_set[sub_face_b])

    # Check whether the two sub_faces are identical in terms of their induced edge crossings
    crossings_a, crossings_b = sub_face_crossings[sub_face_a], sub_face_crossings[sub_face_b]

    # If the sub_faces don't share an edge or aren't identical in edge crossings, they cannot be merged
    if (crossings_a == crossings_b) or (len(merge_edges) == 0):
        return False

    # Determine the new cell's vertex set and edge list
    new_edge_set = sub_face_edge_set[sub_face_a].union(sub_face_edge_set[sub_face_a]) - merge_edges
    new_sub_face = sub_face_a.union(sub_face_b)

    # Update the Cell List
    [sub_faces.remove(sub_face) for sub_face in [sub_face_a, sub_face_b]]
    sub_faces.append(new_sub_face)

    # Update Sub-Face edge list, i.e. remove the two merged ones with the new one
    [sub_face_edge_set.pop(sub_face) for sub_face in [sub_face_a, sub_face_b]]
    sub_face_edge_set[new_sub_face] = new_edge_set

    # Update the Sub-face's induced edge crossing dictionary
    new_sub_face_crossings = copy.deepcopy(sub_face_crossings[sub_face_a])
    [sub_face_crossings.pop(cell) for cell in [sub_face_a, sub_face_b]]
    sub_face_crossings[new_sub_face] = new_sub_face_crossings

    # Update the graph by removing the edges along which the sub-faces were merged
    for merge_edge in merge_edges:
        merge_edge = list(merge_edge)
        graph.remove_edge(u=merge_edge[0], v=merge_edge[1])

    # The merge has successfully occurred
    return True


def select_sub_faces(sub_face_tables, target_faces, target_vertices):

    # Extract sub-faces as dictionary of list; one per face
    sub_faces = {face: sub_face_tables[face].loc[:, "sub_face"].tolist() for face in target_faces}

    # ILP select best pair of sub-faces
    selected_sub_faces, sub_face_incidences = ilp_choose_subface(
        induced_cross_A=sub_face_tables[target_faces[0]].loc[:, target_vertices].to_numpy(),
        induced_cross_B=sub_face_tables[target_faces[1]].loc[:, target_vertices].to_numpy()
    )

    # Convert indices of ILP selection to named faces and vertices
    sub_face_names = tuple([sub_faces[target_faces[i]][selected_sub_faces[i]] for i in range(0, 2)])
    sub_face_neighbors = [[target_vertices[i] for i in range(0, len(target_vertices))
                           if sub_face_incidences[j][i] == 1] for j in range(0, 2)]

    # Store Selection as dictionary
    sub_face_selection = {sub_face_names[i]: frozenset(sub_face_neighbors[i]) for i in range(0, 2)}

    # Return dictionary of selection and their incidence
    return sub_face_selection


def ilp_choose_subface(induced_cross_A, induced_cross_B):

    #
    (W_a, T) = induced_cross_A.shape
    (W_b, T) = induced_cross_B.shape

    # Create model
    m = gp.Model("subfacechoice")

    # Create variables
    # c_i = 1   :   worker i is chosen
    cell_a = m.addVars(W_a, vtype=GRB.BINARY, name="ca")
    cell_b = m.addVars(W_b, vtype=GRB.BINARY, name="cb")

    # e_ij = 1   :   task j is assigned to worker i
    edge_a = m.addVars(W_a, T, vtype=GRB.BINARY, name="ea")
    edge_b = m.addVars(W_b, T, vtype=GRB.BINARY, name="eb")

    # Set objective
    obj = gp.LinExpr(0)
    for j in range(T):
        for i in range(W_a):
            obj.addTerms(induced_cross_A[i][j], edge_a[i, j])
        for i in range(W_b):
            obj.addTerms(induced_cross_B[i][j], edge_b[i, j])
    m.setObjective(obj, GRB.MINIMIZE)

    # Create constraints
    m.addConstr(gp.quicksum(cell_a) <= 1)
    m.addConstr(gp.quicksum(cell_b) <= 1)

    # a task can only be done by one worker
    for j in range(T):
        assignment_amount = gp.LinExpr(0)

        for i in range(W_a):
            assignment_amount.addTerms(1, edge_a[i, j])
        for i in range(W_b):
            assignment_amount.addTerms(1, edge_b[i, j])

        m.addConstr(assignment_amount == 1)

    # a task cannot be assigned to a unchosen worker
    for j in range(T):
        for i in range(W_a):
            m.addConstr(cell_a[i] >= edge_a[i, j])
        for i in range(W_b):
            m.addConstr(cell_b[i] >= edge_b[i, j])

    # Optimize model
    m.optimize()

    index_A, index_B = 0, 0
    subface_A, subface_B = -1, -1

    for v in m.getVars():
        if (v.x) == 1:
            print('%s %g' % (v.varName, v.x))
        if v.varName[0] == "c":
            if v.varName[1] == "a":
                if v.x == 1:
                    subface_A = index_A
                index_A += 1

            if v.varName[1] == "b":
                if v.x == 1:
                    subface_B = index_B
                index_B += 1

    # Get Assignments
    assignment_a = [None] * induced_cross_A.shape[1]
    assignment_b = [None] * induced_cross_B.shape[1]
    regx_numbers = re.compile(r"[+-]?\d+(?:\.\d+)?")
    for v in m.getVars():
        if v.varName[0] != "e":
            continue
        row, column = [int(i) for i in regx_numbers.findall(v.varName)]
        if v.varName[1] == "a" and row == subface_A:
            assignment_a[column] = int(v.x)
        elif v.varName[1] == "b" and row == subface_B:
            assignment_b[column] = int(v.x)

    return (subface_A, subface_B), (assignment_a, assignment_b)


def ilp_choose_face(visibility_matrix):
    (W, T) = visibility_matrix.shape

    # Create model
    m = gp.Model("facechoice")

    # Create variables
    # c_i = 1   :   worker i is chosen
    cell = m.addVars(W, vtype=GRB.BINARY, name="c")

    # e_ij = 1   :   task j is assigned to worker i
    edge = m.addVars(W, T, vtype=GRB.BINARY, name="e")

    # Set objective
    obj = gp.quicksum(edge)
    m.setObjective(obj, GRB.MAXIMIZE)

    # Create constraints
    # can only select two workers
    m.addConstr(gp.quicksum(cell) <= 2)

    # a task can only be done by one worker
    for j in range(T):
        assignment_amount = gp.LinExpr(0)

        for i in range(W):
            assignment_amount.addTerms(1, edge[i, j])

        m.addConstr(assignment_amount <= 1)

    # a task cannot be assigned to a unchosen worker
    for i in range(W):
        for j in range(T):
            m.addConstr(cell[i] >= edge[i, j])

    # if a worker cant do the task, it is not assigned to him
    for i in range(W):
        for j in range(T):
            m.addConstr(edge[i, j] <= visibility_matrix[i][j])

    # Optimize model
    m.optimize()

    faces = []
    index = 0

    for v in m.getVars():

        if v.varName[0] == "c":

            if (v.x) == 1:
                faces += [index]
        index += 1

    return (faces)


def find_all_subfaces(target_faces, face_vertex_map, graph):

    # TODO: rework this nonsense to work with our updated data structures and
    # Prepare dictionary of sets within which to store all found faces per face
    subfaces = {face: None for face in target_faces}

    # Limit search of cycle basis to subfaces
    for face in target_faces:

        # Extrac
        vertex_set = face_vertex_map.get(face)

        # Keep only vertices pertaining to current subgraph and search for cycle basis within
        subgraph = copy.deepcopy(graph).subgraph(list(vertex_set))
        subfaces[face] = find_inner_faces(subgraph)  # todo: also pass positions?

    # Return all subfaces for each target face
    return subfaces


def get_ordered_face_edges(faces, graph):
    ordered_face_edges = dict.fromkeys(faces)
    for face in faces:
        vertex_sequence = list(face) if len(face) == 1 else get_face_vertex_sequence(face, graph)
        ordered_face_edges[face] = vertex_sequence
    return ordered_face_edges


def unlist(nested_list):
    return list(it.chain.from_iterable(nested_list))


def update_faces_with_edge_map(face_incidence_table, face_edge_map, edge_map):

    for index, row in face_incidence_table.iterrows():
        face = row["identifier"]
        face_edges = face_edge_map[face]
        new_face_edges = []
        for edge in face_edges:
            new_face_edges.append(edge_map.get(edge, [edge]))
        new_face_edges = unlist(new_face_edges)
        new_face_identifier = frozenset(unlist(new_face_edges))

        face_edge_map.pop(face)
        face_edge_map[new_face_identifier] = new_face_edges
        face_incidence_table.at[index, "identifier"] = new_face_identifier


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
                print("SHRINKING CYCLES are fucked")
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
            sys.exit("CYCLES FUCKED")

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


def find_inner_faces(graph, positions=None, as_set=True):

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
    cycle_ordered_edges = get_ordered_face_edges(cycles, graph)

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


def find_outer_face(ordered_face_edges, graph, positions):
    # TODO: identify "faces" which consist of non-cycles, i.e. disconnected "lines" of vertices

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
                # TODO: see above todo!
                sys.exit(f"Edge ({edge_a}, {edge_b}) did not map to a face.")

    # Find edges which only map to a single face
    edges = set([frozenset(edge) for edge in faces_per_edge.keys() if faces_per_edge[edge] == 1])

    # Identify Unique sets of edges to find faces
    faces = copy.deepcopy(edges)
    find_vertex_sets_from_edges(faces)

    # Map whether face is a cycle or not
    face_is_cycle = {face: True for face in faces}

    # Map edges to faces
    face_edge_sets = find_face_edge_sets(faces, edges)

    # Find singleton vertices and Check whether singleton falls outside for
    vertices = [vertex for vertex in graph.nodes if len(graph.edges(vertex)) == 0]
    for face in faces:
        sorted_vertices = get_sorted_face_vertices([tuple(edge) for edge in face_edge_sets[face]], is_sorted=False)
        cycle_coordinates = [positions[vertex] for vertex in sorted_vertices]
        cycle_path = mpltPath.Path(cycle_coordinates[0:-1])
        [vertices.remove(v) for v in copy.copy(vertices) if cycle_path.contains_points([positions[v]])]

    # Update Edges Sets and Face with Singletons in the outer face
    vertices = set([frozenset([vertex]) for vertex in vertices])
    faces = faces.union(vertices)
    [face_is_cycle.update({vertex: False}) for vertex in vertices]

    print(f"outer faces: {faces}")
    # Return unique vertex sets from the found singleton edges
    return faces, face_edge_sets, face_is_cycle


def find_face_edge_sets(faces, edges):
    face_edge_sets = {face: set() for face in faces}
    [face_edge_sets[face].add(edge) for edge in edges for face in faces if len(edge.intersection(face)) == 2]
    return face_edge_sets


def find_vertex_sets_from_edges(edge_sets):
    for edge_set_a in edge_sets:
        for edge_set_b in edge_sets:
            if edge_set_a == edge_set_b:
                continue
            if not edge_set_a.intersection(edge_set_b):
                continue

            new_edge_set = edge_set_a.union(edge_set_b)
            edge_sets.remove(edge_set_a)
            edge_sets.remove(edge_set_b)
            edge_sets.add(new_edge_set)

            return find_vertex_sets_from_edges(edge_sets)


def build_face_to_edge_map(graph, faces):

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
    # Initialize target vertices as set
    target_vertex_set = set(target_vertices)

    # Initialize an empty dictionary (using sets as keys) to store vertex incidence sets
    face_incidences = {face: face & target_vertex_set for face in faces}

    # Return Dictionary of face vertex incidence
    return face_incidences


def find_outer_face_vertex_incidence(outer_face, inner_faces, target_vertices):

    # Initialize an empty dictionary (using sets as keys) to store vertex incidence sets
    outer_face_incidences = {outer_face: {}}  # face_incidences = {face: dict() for face in faces}
    outer_face_incidence = outer_face.intersection(target_vertices)

    # Initialize set of vertex set as list, and target vertices as set
    target_vertex_set = set(target_vertices)

    # Iterate over all faces
    for inner_face in inner_faces:

        # Determine how many target vertices are incident and left over
        inner_face_incidence = inner_face.intersection(target_vertex_set)
        outer_face_incidences[outer_face][inner_face] = outer_face_incidence.union(inner_face_incidence)

    # Return Dictionary of face vertex incidence
    return outer_face_incidences


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


def is_inner_face_convex(ordered_face_edges, positions):
    face_vertices = get_sorted_face_vertices(ordered_face_edges, is_sorted=True)
    face_angles = calculate_face_inner_angles(face_vertices, positions)
    return all([angle <= 180 for angle in face_angles])
