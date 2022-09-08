

def get_face_vertex_sequence(face, graph):
    # todo: does not work if edges exist between elements of the face other than the minimal cycle
    face_edges = [None] * len(face)
    current = 0
    for edge in graph.edges:
        common_vertices = face.intersection(set(edge))
        if len(common_vertices) == 2:
            face_edges[current] = edge
            current += 1
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
