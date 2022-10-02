from src.edge_crossings import *
from src.edges import *

import itertools as it
import networkx as nx
import numpy as np
import copy


def get_nested_key(dictionary, key_a, key_b):
    try:
        return dictionary[key_a][key_b]
    except KeyError:
        return None


def draw_all_line_segments(graph, positions, virtual_edge_set, bounds, already_extended, start_index=None):

    # Create new graph and positions objects
    segment_graph, segment_positions = copy.deepcopy(graph), copy.deepcopy(positions)
    segment_edge_map = dict()

    # Store nodes and edges for easier look-up
    # virtual_nodes, corner_nodes = nx.get_node_attributes(graph, "virtual"), nx.get_node_attributes(graph, "corner")
    # real_nodes = [v for v in segment_graph.nodes if (virtual_nodes.get(v, 0) != 1) and (corner_nodes.get(v, 0) != 1)]
    real_nodes = [node for node, real in nx.get_node_attributes(graph, "real").items() if real == 1]
    print(f"real nodes: {real_nodes}")
    edges = frozenset([frozenset(edge) for edge in list(segment_graph.edges())])

    # Store the number of nodes and largest index
    vertex_index = start_index if start_index else max(graph.nodes)
    number_of_nodes = len(real_nodes)

    # Iterate over all pairwise vertex combinations
    for index_a in range(0, number_of_nodes):
        vertex_a = real_nodes[index_a]
        for index_b in range(index_a + 1, number_of_nodes):
            vertex_b = real_nodes[index_b]

            # Calculate intersections with boundary
            intersections = extend_line(segment_positions[vertex_a], segment_positions[vertex_b], bounds)

            # # Check if this particular combination of vertices is connected (by virtual edge sets)
            virtual_connection = [v_edge_set for v_edge_set in virtual_edge_set if {vertex_a, vertex_b} <= v_edge_set]
            if virtual_connection:
                virtual_edges = virtual_edge_set[virtual_connection[0]]
                [segment_graph.add_edge(u_of_edge=e[0], v_of_edge=e[1], segment=1) for e in virtual_edges
                 if not segment_graph.has_edge(e[0], e[1])]
                already_connected = 1
            else:
                already_connected = 1 if {real_nodes[index_a], real_nodes[index_b]} in edges else 0
                virtual_edges = [{real_nodes[index_a], real_nodes[index_b]}] if already_connected else []

            # If the two vertices are not yet connected, connect them with a segment edge
            if already_connected == 0:
                segment_graph.add_edge(u_of_edge=vertex_a, v_of_edge=vertex_b, segment=1)
                virtual_edges = [{vertex_a, vertex_b}]

            # Project vertex-a onto vertex-b and vice-versa
            new_boundary_vertices = [None, None]
            for pair_index, (target, joint) in enumerate([(vertex_b, vertex_a), (vertex_a, vertex_b)]):

                # Check if vertex pair has already been extended
                extended_vertex = get_nested_key(already_extended, target, joint)
                if extended_vertex is not None:

                    # Check if any virtual edges exist between joint and extension terminus
                    extended_edges = virtual_edge_set.get(frozenset((joint, extended_vertex[0])), None)
                    if extended_edges is not None:
                        print(f"\nextended edges for {vertex_a} and {vertex_b}: {extended_edges}")
                        print(f"joint: {joint}")
                        connections = get_vertex_sequence(edges=extended_edges,
                                                          first_node=joint,
                                                          is_ordered=False)

                    # If no additional vertices exist, then connect only the joint and the
                    else:
                        connections = [joint, extended_vertex[0]]

                # If no connections exist, the keep and empty list
                else:
                    connections = [joint]

                # Add another virtual vertex beyond the boundary
                vertex_index += 1
                new_boundary_vertices[pair_index] = vertex_index
                segment_graph.add_node(node_for_adding=vertex_index, boundary=1)
                segment_positions[vertex_index] = np.asarray(intersections[pair_index])
                connections.append(vertex_index)

                # Add new edges (line segments) to graph if the edge does not yet exist
                print(f"connections: {connections}")
                for index in range(1, len(connections)):
                    # edge = {connections[index - 1], connections[index]}
                    connection_a, connection_b = connections[index - 1], connections[index]
                    virtual_edges.append((connection_a, connection_b))
                    if {connection_a, connection_b} in edges:
                        print(f"{connection_a, connection_b} already in edges")
                        continue
                    segment_graph.add_edge(u_of_edge=connection_a, v_of_edge=connection_b, segment=1)

            # Append to new edge map
            print(f"virtual_edges: {virtual_edges}")
            segment_edge_map[frozenset(new_boundary_vertices)] = [set(edge) for edge in virtual_edges]

    print(f"\nSEGMENT MAP:")
    [print(f"{key} - {item}") for key, item in segment_edge_map.items()],

    # Return new graph and positions objects
    return segment_graph, segment_positions, segment_edge_map


def do_both_vertices_map_to_either(vertex_a, vertex_b, dict_a, dict_b):
    return all([does_one_vertex_map(vertex_a, vertex_b, dictionary) for dictionary in [dict_a, dict_b]])


def does_one_vertex_map(vertex_a, vertex_b, dictionary):
    return any([dictionary.get(v, 0) == 1 for v in (vertex_a, vertex_b)])


def do_both_vertices_map(vertex_a, vertex_b, dictionary):
    return all([dictionary.get(v, 0) == 1 for v in (vertex_a, vertex_b)])


def cull_by_vertex_identity(vertex_a, vertex_b, face):
    return len(face.intersection({vertex_a, vertex_b})) == 2


def cull_all_line_segment_graph(target_faces, face_edge_map, face_vertex_map, segment_edge_map, graph, positions):

    # Create new graph objects
    culled_graph, culled_positions = copy.deepcopy(graph), copy.deepcopy(positions)
    virtual_edges = [set(e) for e, v in nx.get_edge_attributes(G=culled_graph, name="virtual").items() if v == 1]
    print(f"virtual edges: {virtual_edges}")

    # Keep track of non-intersecting segments
    segment_edges = [e for e, v in nx.get_edge_attributes(G=culled_graph, name="segment").items() if v == 1]
    edges_to_keep = set()

    # Initialize empty nested dictionary
    face_intersection_map = {target_face: dict() for target_face in target_faces}

    #
    for face in target_faces:

        #
        for segment_edge, sub_edges in segment_edge_map.items():

            #
            found_in_vertex_set = False

            # Try mapping
            for (vertex_a, vertex_b) in sub_edges:

                if len(face_vertex_map[face].intersection({vertex_a, vertex_b})) == 2:
                    edges_to_keep.add(frozenset((vertex_a, vertex_b)))
                    found_in_vertex_set = True

            #
            if found_in_vertex_set:
                continue

            #
            for (vertex_a, vertex_b) in sub_edges:

                intersections_found = 0

                # Convert edge to frozenset
                edge = frozenset((vertex_a, vertex_b))

                # Skip if the both edge vertices are virtual
                if edge in virtual_edges:
                    continue

                # Check Whether Edge Crosses the faces
                for face_edge in face_edge_map[face]:

                    # Determine if line actually passes through one of the vertices which define the face
                    if {vertex_a, vertex_b} & {face_edge[0], face_edge[1]}:
                        intersection_key = face_edge[0] if face_edge[0] in {vertex_a, vertex_b} else face_edge[1]
                        intersection = culled_positions[intersection_key]

                    # Check edge crossings
                    else:
                        intersection_key = (face_edge[0], face_edge[1])
                        edge_point_a, edge_point_b = culled_positions[vertex_a], culled_positions[vertex_b]
                        face_point_a, face_point_b = culled_positions[face_edge[0]], culled_positions[face_edge[1]]
                        intersection = line_intersection(edge_point_a, edge_point_b, face_point_a, face_point_b)

                    # Store intersection if one exists
                    if intersection is not None:
                        if edge not in face_intersection_map[face].keys():
                            face_intersection_map[face][edge] = dict()
                        if intersection_key not in face_intersection_map[face][edge].keys():
                            intersections_found += 1
                            face_intersection_map[face][edge][intersection_key] = intersection

                # If only one intersection was found, remove the edge -> not sufficient numbers
                if intersections_found == 1:
                    face_intersection_map[face].pop(edge)

    # Remove Edges and Vertices which did not intersect any face
    for edge in segment_edges:
        intersected = [frozenset(edge) in face_intersection_map[face].keys() for face in target_faces]
        if frozenset(edge) in edges_to_keep or any(intersected):
            continue
        culled_graph.remove_edge(u=edge[0], v=edge[1])

    # Remove Edges that are not segments, but sight cell defining edges; first track real and target sight cell edges
    real_edges = [set(e) for e, v in nx.get_edge_attributes(G=culled_graph, name="real").items() if v == 1]
    target_edges = set([frozenset(e) for target_face in target_faces for e in face_edge_map.get(target_face)])

    # Iterate over all sight cells that are not target sight cells
    for face, face_edges in face_edge_map.items():
        if face in target_faces:
            continue

        # Iterate over all edges in these sight cells
        for face_edge in face_edges:

            # Edge may already have been removed, if so, skip
            if not culled_graph.has_edge(u=face_edge[0], v=face_edge[1]):
                continue

            # Edge may be real, i.e. bound defining but part of original graph
            if {face_edge[0], face_edge[1]} in real_edges:
                continue

            # Edge may also be shared by target sight cells
            if {face_edge[0], face_edge[1]} in target_edges:
                continue

            # If none of the above conditions are met, remove the edge
            culled_graph.remove_edge(u=face_edge[0], v=face_edge[1])

    # Remove all vertices with degree zero
    no_degree_nodes = [v for v in culled_graph.nodes if culled_graph.degree(v) == 0]
    [culled_graph.remove_node(v) for v in no_degree_nodes]

    # Return new graph, positions, and intersection map
    return culled_graph, culled_positions, face_intersection_map


def create_subface_graph(graph, positions, target_faces, face_vertex_map, face_intersection_map):

    # Iterate over all (hopefully 2) target faces
    node_list = list(graph.nodes())
    vertex_index = max(node_list)

    edges_to_be_removed = set()
    nodes_to_be_removed = [n for n, b in nx.get_node_attributes(G=graph, name="boundary").items() if b == 1]
    edge_to_virtual_vertex = dict()
    # face_vertex_map = {face: set() for face in target_faces}

    # Update face vertex map to feature sets, not frozensets
    face_vertex_map.update({face: set(vertices) for face, vertices in face_vertex_map.items()})

    for target_face in target_faces:
        # [face_vertex_map[target_face].add(face_vertex) for face_vertex in list(target_face)]

        for intersecting_edge in face_intersection_map[target_face].keys():
            edges_to_be_removed.add(frozenset(intersecting_edge))
            edge_targets = []
            for face_edge in face_intersection_map[target_face][intersecting_edge].keys():

                # Ensure the intersection is with a line, not a face vertex
                if type(face_edge) is tuple:

                    # Add currently considered face edge to the dictionary
                    if face_edge not in edge_to_virtual_vertex:
                        edge_to_virtual_vertex[face_edge] = set()
                    edges_to_be_removed.add(frozenset(face_edge))
                    intersection = face_intersection_map[target_face][intersecting_edge][face_edge]

                    # Create new Vertex
                    vertex_index += 1
                    edge_targets.append(vertex_index)
                    graph.add_node(node_for_adding=vertex_index, segment=1)
                    face_vertex_map[target_face].add(vertex_index)
                    positions[vertex_index] = np.asarray(intersection)
                    edge_to_virtual_vertex[face_edge].add(vertex_index)
                else:
                    edge_targets.append(face_edge)

            #
            graph.add_edge(u_of_edge=edge_targets[0], v_of_edge=edge_targets[1], segment=1)

    # Add virtual edge connections
    # [print(f"{index} - {edge_to_virtual_vertex[index]}") for index in edge_to_virtual_vertex.keys()]
    virtual_edge_set = add_virtual_edges(graph, positions, edge_to_virtual_vertex)

    #
    subface_vertex_map = {face: face_vertex_map[face] for face in target_faces}

    # Remove
    [graph.remove_edge(u=list(edge)[0], v=list(edge)[1]) for edge in edges_to_be_removed]
    [graph.remove_node(vertex) for vertex in nodes_to_be_removed]

    return virtual_edge_set, subface_vertex_map


def update_face_vertex_map(vertex_map, virtual_edge_map):

    # Lambda function to turn a nested list into a single list
    unlist = lambda nested_list: list(it.chain.from_iterable(nested_list))

    # Grow list of all vertices which are featured in virtual edges
    new_vertices = {face: [] for face in vertex_map.keys()}

    # Iterate over all faces
    for face, vertices in vertex_map.items():

        # Iterate over all edges and virtual edges from the planarization
        for edge, virtual_edges in virtual_edge_map.items():

            # If no virtual edges exist for an edge, skip it
            if not virtual_edges:
                continue

            # Check if both vertices of the edge are part of the current face's vertex set
            if len(vertices.intersection(set(edge))) == 2:

                # Grow the list of vertices featured in virtual edges
                new_vertices[face] += unlist(virtual_edges)

        # If no edge crossings exist and hence no new vertices, move on to the next face
        if not new_vertices[face]:
            continue

        # Update the face's vertex set as the union of the old set and the new set of edge crossing vertices
        vertex_map[face] = vertices.union(set(new_vertices[face]))
