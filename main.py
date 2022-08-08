from src.graph_simulation import *
from src.graph_drawing import *
from src.edge_crossings import *
from src.faces import *

import numpy as np

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Define Input Parameters
    embedding = "kamada_kawai"
    n_vertices = 20
    m_edges = 3
    seed = 1

    # EMBED INPUT GRAPH ------------------------------------------------------------------------------------------------

    # Create Simulated Graph
    graph = create_barabasi_albert_graph(n=n_vertices, m=m_edges, seed=seed)
    positions = embed_graph(graph=graph, embedding=embedding)
    virtual_index_start = graph.number_of_nodes()

    # Collect and Check Edge Crossings
    edge_crossings, vertex_crossings = locate_edge_crossings(graph, positions)
    assert vertex_edge_crossing_equality(vertex_crossings, edge_crossings), \
        "Vertex and Edge Crossing Numbers not equivalent. Sum(Cr(Vi)) / 4 = Sum(Cr(Ej))"

    # Draw and Save Non-Planar rGraph
    draw_graph(graph=graph, positions=positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=0)
    save_drawn_graph(output_path)

    for vertex in positions.keys():
        print("{} - {} : {}".format(vertex, graph.nodes[vertex], positions[vertex]))

    # LOCATE AND REMOVE TARGET VERTEX ----------------------------------------------------------------------------------

    # Find vertex involved in the largest number of edge crossings
    target_vertex_index = get_target_vertex_index(vertex_crossings, graph)
    target_vertex_adjacency = list(graph.neighbors(target_vertex_index))
    print("Split Target Vertex = Vertex #{}".format(target_vertex_index))
    print("Target Adjacency    = {}".format(target_vertex_adjacency))

    # Delete Target Vertex from Graph
    print(edge_crossings)
    remaining_edge_crossings = get_remaining_edge_crossings(graph, edge_crossings, target_vertex_index)
    remaining_graph, remaining_positions = remove_target_vertex(graph, positions, target_vertex_index)

    # Draw and Save Graph with target removed
    draw_graph(graph=remaining_graph, positions=remaining_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=1)
    save_drawn_graph(output_path)

    for vertex in remaining_positions.keys():
        print("{} - {} : {}".format(vertex, remaining_graph.nodes[vertex], remaining_positions[vertex]))

    # DEBUG ------------------------------------------------------------------------------------------------------------

    # Add Nodes for each edge crossing
    debug_graph = copy.deepcopy(remaining_graph)
    debug_positions = copy.deepcopy(remaining_positions)
    n = virtual_index_start
    for edge_a in remaining_edge_crossings.keys():
        for edge_b in remaining_edge_crossings[edge_a].keys():
            debug_graph.add_node(node_for_adding=n, split=0, target=1, virtual=0)
            debug_positions[n] = np.asarray(remaining_edge_crossings[edge_a][edge_b])
            n += 1

    # Draw and Save Graph with target removed
    draw_graph(graph=debug_graph, positions=debug_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=2)
    save_drawn_graph(output_path)

    # PLANARIZATION ----------------------------------------------------------------------------------------------------

    # Planarize Graph
    plane_graph, plane_positions = planarize_graph(graph=remaining_graph,
                                                   positions=remaining_positions,
                                                   edge_crossings=remaining_edge_crossings,
                                                   starting_index=virtual_index_start)
    planar_edge_crossings, planar_vertex_crossings = locate_edge_crossings(plane_graph, plane_positions)

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=3)
    save_drawn_graph(output_path)

    # FACE IDENTIFICATION ----------------------------------------------------------------------------------------------

    # Locate faces and best two for target face
    faces = frozenset(find_all_faces(plane_graph))
    face_edge_map = build_face_to_edge_map(plane_graph, faces)
    # TODO: FIX THE MAPPING TO NOT INCLUDE SELF-LOOPS (SEE Picture 4)

    face_incidences = find_face_vertex_incidence(faces, target_vertex_adjacency)
    max_incidence, selected_faces = get_maximally_incident_faces(face_incidences)
    print("#Face Target Sets Found: {}".format(len(selected_faces)))
    print("Selected Faces: {}".format(selected_faces[1]))
    # TODO: RE-EVALUATE HOW LEGAL THE CREATED AND SELECTED FACES ARE

    for face in selected_faces[1]:
        print("Face: {}".format(face))
        for edge in face_edge_map[frozenset(face)]:
            print("Edge: {}".format(edge))
            plane_graph.edges[edge]["target"] = 1
    print(face_incidences[selected_faces[0][0]][selected_faces[0][1]])
    print(selected_faces)

    # Draw and Save Planar rGraph
    draw_graph(graph=plane_graph, positions=plane_positions)
    output_path = create_output_path(embedding=embedding, n_vertices=n_vertices, m_edges=m_edges, seed=seed, n_splits=4)
    save_drawn_graph(output_path)


