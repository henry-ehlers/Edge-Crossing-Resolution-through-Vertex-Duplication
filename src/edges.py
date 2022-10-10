import itertools as it
import copy


def unlist(nested_list):
    return list(it.chain.from_iterable(nested_list))


def get_valid_starting_vertex(edges: [(int, int)], starting_vertex: int = None) -> int:
    """
    A function to either obtain a valid starting vertex from a list of edges, or ensure the provided starting vertex is
    indeed a valid choice, i.e. that the provided list of edges is either a cycle of edges or a sequence with two
    defined end-points
    :param edges: a list of (unsorted) edges as tuples to be sorted
    :param starting_vertex: an optional single integer vertex at which to start the sorting
    :return: the integer identifier of the starting vertex
    """

    # Extract a list of vertices from the list of edges
    vertices = unlist(edges)

    # If a starting vertex is provided, ensure it actually exists in the vertex list
    if starting_vertex is not None:
        assert starting_vertex in vertices, \
            f"Providing starting vertex {starting_vertex} is not in edge's vertex set {vertices}."

    # Count the number of occurences of each vertex in the list of vertices
    vertex_counts = {vertex: 0 for vertex in list(set(vertices))}
    for vertex in vertices:
        vertex_counts[vertex] += 1

    # Count the number of "knots" and "forks" in the edge sequence
    forks = [node for node, count in vertex_counts.items() if count > 2]
    assert len(forks) == 0, \
        f"Edge Sequence has no unique order, as nodes {forks} are featured more than twice."

    # Count the number of ends in the sequence is only 2, or all vertices must appear exactly twice
    ends = [node for node, count in vertex_counts.items() if count == 1]
    cycles = [node for node, count in vertex_counts.items() if count == 2]
    assert (len(ends) == 2) or (len(cycles) == len(vertex_counts)), \
        f"Edge sequence has too many or too few end vertices, as nodes {ends} all only appear once, or the number" \
        f"of double matched vertices, {len(cycles)} is not the length of the available nodes {len(vertex_counts)}"

    # If the starting vertex is provided, ensure that it corresponds to end point vertex
    if (starting_vertex is not None) and (len(ends) > 0):
        assert starting_vertex in ends, \
            f"Provided starting vertex {starting_vertex} not a valid selection from end points {ends}."

    # Return a starting vertex if it was provided and is legal
    if starting_vertex is not None:
        return starting_vertex

    # If a legal sequence of edges was provided, return one of its end-points
    elif len(ends) > 0:
        return ends[0]

    # If a legal cycle was provided, return an (effectively) random vertex
    else:
        return vertices[0]


def get_first_edge(edges: [(int, int)], starting_vertex: int) -> (int, int):
    """
    A function which returns an edge which features the provided starting vertex in its first position
    :param edges: a list of (unsorted) edges as tuples to be sorted
    :param starting_vertex: an optional single integer vertex at which to start the sorting
    :return: a sorted list of edges
    """

    # Iterate over all edges in set
    for edge in edges:

        # Check whether first vertex of edge is the starting vertex
        if edge[0] == starting_vertex:
            return edge

        # Check whether second vertex of edge is the starting vertex
        elif edge[1] == starting_vertex:
            return edge[1], edge[0]

    # Provided starting vertex Mapped to no edge
    assert True, f"Starting Vertex {starting_vertex} mapped to no edge {edges}"


def get_ordered_edges(edges, starting_vertex: int = None) -> [(int, int)]:
    """
    A function to sort a list of edges provided as tuples, such that the second vertex of an edge corresponds to the
    first vertex of the next edge, i.e. (a, c) (c, b) (b, a)
    :param edges: a list of (unsorted) edges as tuples to be sorted
    :param starting_vertex: an optional single integer vertex at which to start the sorting
    :return: a sorted list of edges
    """

    # Convert input edges to tuples to accommodate set types
    edges = [edge if type(edge) is tuple else tuple(edge) for edge in edges]
    print(f"edges: {edges}")
    # Get a valid starting vertex / Ensure the provided one is valid
    starting_vertex = get_valid_starting_vertex(edges, starting_vertex)

    # Convert list of tuples to list of frozen-sets
    remaining_edges = [frozenset(edge) for edge in edges]

    # Obtain an edge which features the starting vertex
    first_edge = get_first_edge(edges, starting_vertex)

    # Initialize the list of sorted edges, and remove first one from the remaining set
    sorted_edges = [first_edge] + [(None, None)] * (len(edges) - 1)
    remaining_edges.remove(frozenset(first_edge))

    # Iterate over each remaining index to be sorted
    for i in range(1, len(sorted_edges)):

        # Iterate over all remaining edges to be included
        for edge_set in remaining_edges:

            # Convert Set Edge to forward and 'reversed' edge
            f_edge = tuple(edge_set)
            r_edge = f_edge[::-1]

            # Check whether the forward facing edge matches the last sorted edge
            if f_edge[0] == sorted_edges[i - 1][1]:
                sorted_edges[i] = f_edge
                remaining_edges.remove(frozenset(f_edge))
                break

            # Check whether the 'reversed' facing edge matches the last sorted edge
            elif r_edge[0] == sorted_edges[i - 1][1]:
                sorted_edges[i] = r_edge
                remaining_edges.remove(frozenset(r_edge))
                break

    # Ensure no edges remain to be sorted
    assert len(remaining_edges) == 0, \
        f"Not all edges sorted. These remain: {sorted_edges}"

    # Return sorted edges as list of tuples
    return sorted_edges


def get_vertex_sequence(edges, starting_vertex=None, is_ordered=False):

    # If the edges are ordered, they may not need to be sorted
    if is_ordered:

        # If no starting vertex is provided, the ordered edges can be taken as is
        if starting_vertex is None:
            ordered_edges = copy.deepcopy(edges)

        # If a starting vertex is provided, ensure it corresponds to the first vertex of the edges
        else:

            # If it does, the ordered edges can be utilized as is
            if edges[0][0] == starting_vertex:
                ordered_edges = copy.deepcopy(edges)

            # If not, re-sort the provided edges with the provided starting vertex
            else:
                ordered_edges = get_ordered_edges(edges=edges, starting_vertex=starting_vertex)

    # If the provided edges are not sorted, sort them using the provided starting vertex
    else:
        ordered_edges = get_ordered_edges(edges=edges, starting_vertex=starting_vertex)

    # Check whether the ordered edge list's first and last vertices are the same -> it's a cycle
    if ordered_edges[0][0] == ordered_edges[-1][1]:
        vertex_sequence = [edge[0] for edge in ordered_edges]

    # If it is not a cycle, the end vertices must be treated differently
    else:
        vertex_sequence = [edge[0] for edge in ordered_edges] + ordered_edges[-1][1]

    # Return the ordered vertex sequence
    return vertex_sequence


def get_face_vertex_sequence(face, graph):
    # todo: does not work if edges exist between elements of the face other than the minimal cycle
    #face_edges = [None] * len(face)
    face_edges = []
    # print(f"face: {face}")
    for edge in graph.edges:
        # print(f"edge: {edge}")
        common_vertices = face.intersection(set(edge))
        if len(common_vertices) == 2:  # TODO: this does not always produce the edge sequence if there is a triangle
            face_edges.append(edge)
    # print(f"face edges: {face_edges}")
    sorted_face_edges = sort_face_edges(face_edges)
    return sorted_face_edges


def sort_face_edges(edge_list):

    # TODO: reformulate search in terms of indeces that can be 'blacked' out because there were already included

    # Convert set of frozensets to list of tuples
    # TODO: make this check more robust
    edge_list = [tuple(edge) for edge in edge_list] if isinstance(edge_list, set) else edge_list

    # Initialize new list of sorted edges in a cycle
    new_list = [(None, None)] * len(edge_list)
    new_list[0] = edge_list[0]

    # Iterate over all indices to be filled
    for index in range(1, len(edge_list)):

        # Specify the next target value as the last vertex of the first edge
        target_value = new_list[index-1][1]

        # Find the next possible edge which matches the last edge's second vertex and has not yet been included
        next_element = [element for element in edge_list if (target_value in element) and
                        (element not in new_list) and ((element[1], element[0]) not in new_list)][0]

        # Store either the original edge or its reverse depending on which allows for the cycle to continue
        reverse_edge = (next_element[1], next_element[0])
        new_list[index] = reverse_edge if reverse_edge[0] == target_value else next_element

    # Ensure that the found order is indeed cyclical
    assert new_list[0][0] == new_list[-1][1], \
        "Found sorted order of edges is not cyclical."

    # Return the sorted cycle
    return new_list


def get_sorted_face_vertices(edges, is_sorted=False):
    if not is_sorted:
        if not all(edges[index][1] == edges[(index + 1) % len(edges)] for index in range(0, len(edges))):
            edges = sort_face_edges(edges)
    return [edge[0] for edge in edges]
