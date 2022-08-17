import networkx as nx


def get_face_vertex_sequence(face, graph):
    face_edges = [None] * len(face)
    current = 0
    for edge in graph.edges:
        if len(set(edge) & face) == 2:
            face_edges[current] = edge
            current += 1
    sorted_face_edges = sort_edges(face_edges)
    return sorted_face_edges


def sort_edges(edge_list):
    input_dict = {edge[0]: edge[1] for edge in edge_list}
    input_dict.update({edge[1]: edge[0] for edge in edge_list})
    elem = edge_list[0][0]  # start point in the new list
    new_list = []  # List of tuples for holding the values in required order
    for _ in range(len(edge_list)):
        possibility_1 = (elem, input_dict[elem])
        possibility_2 = (input_dict[elem], elem)
        new_list.append(possibility_1) if possibility_1 not in new_list else new_list.append(possibility_2)
        elem = input_dict[elem]
    return new_list

