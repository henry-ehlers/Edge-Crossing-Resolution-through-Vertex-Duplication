from src.create_random_graphs import create_barabasi_albert_graph

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    simulated_graph = create_barabasi_albert_graph(n=20, m=5, seed=1)
    print(simulated_graph)

