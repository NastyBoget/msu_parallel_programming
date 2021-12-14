import numpy as np
import os
import seaborn as sns


edges_num = 0


def print_edges(key1: str, key2: str, graph_dict: dict, f, colors) -> None:
    global edges_num
    if key1 not in graph_dict or key2 not in graph_dict:
        print(f"{key1} or {key2} not in graph")
        return
    if graph_dict[key1][2] == graph_dict[key2][2]:
        print(f'"{key1}" -- "{key2}" [color="{graph_dict[key1][2]}",penwidth=0.5]', file=f)
    else:
        edges_num += 1
        print(f'"{key1}" -- "{key2}" [color="{colors[-1]}",penwidth=0.5]', file=f)


def visualize(file_name: str) -> None:
    with open(file_name) as f:
        lines = f.readlines()
    lines = [list(map(float, line.strip().split())) for line in lines]
    lines = np.array(lines)

    n = int(max(lines[:,-1]))
    colors = sns.color_palette(None, n + 2)
    colors = colors.as_hex()

    graph_dict = {}
    for line in lines:
        graph_dict[f"{int(line[0])} {int(line[1])}"] = [line[2], line[3], colors[int(line[-1])]]

    with open(f"{file_name}.dot", "w") as f:
        print("graph {", file=f)
        print("graph [ dpi = 500 ]; ", file=f)
        print("node [shape=point, height=0.01]", file=f)
        for key, value in graph_dict.items():
            print(f'"{key}" [pos="{value[0]},{value[1]}!", color="{value[2]}"]', file=f)
        max_i = int(max(lines[:, 0]))
        max_j = int(max(lines[:, 1]))
        for i in range(max_i + 1):
            for j in range(max_j + 1):
                if i < max_i:
                    key1, key2 = f"{i} {j}", f"{i + 1} {j}"
                    print_edges(key1, key2, graph_dict, f, colors)
                if j < max_j:
                    key1, key2 = f"{i} {j}", f"{i} {j + 1}"
                    print_edges(key1, key2, graph_dict, f, colors)
        print("}", file=f)
    print(f"Number of edges: {edges_num}")
    os.system(f"dot -Kneato -n -Tpng -O {file_name}.dot")
    os.system(f"rm -rf {file_name}.dot")


if __name__ == "__main__":
    visualize("out.txt")
