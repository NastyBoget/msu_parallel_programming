import matplotlib.pyplot as plt
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

    n = int(max(lines[:, -1]))
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


def make_graphics(statistics: dict, kind: int) -> None:
    text_dict = {"t": "Время выполнения", "s": "Ускорение", "e": "Эффективность"}
    color_dict = {"t": '-bD', "s": '-gD', "e": '-rD'}

    for key in ("t", "s", "e"):
        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(212)
        axes = [ax1, ax2, ax3]
        fig.set_figheight(15)
        fig.set_figwidth(15)
        for k, n in enumerate(statistics):
            x = statistics[n]["p"]
            y = np.array(statistics[n][key])

            axes[k].set(title="{}: результаты для сетки {}x{}".format(text_dict[key], int(n), int(n)),
                        xlabel='Число процессов', ylabel=text_dict[key])

            axes[k].plot(x, y, color_dict[key], markevery=np.arange(len(y)))
            for i_x in range(len(x)):
                axes[k].annotate(str((x[i_x], round(y[i_x], 3))), (x[i_x], y[i_x]))

        fig.savefig('images/{}_{}.jpg'.format(kind, key))


# K, NET - fixed
# S = T1 / T
# E = T1 / (T * p)
# ------------------------------------------------
# | num_processes | net_size | T | S | E | edges |
# ================================================
# |       1       | 128x128  |   |   |   |       |
# ------------------------------------------------
# ...
# ------------------------------------------------
# |      256      | 128x128  |   |   |   |       |
# ================================================
# |       1       | 256x256  |   |   |   |       |
# ------------------------------------------------
# ...

def add_statistics(result_dict_n: dict) -> None:
    s_list = []
    e_list = []
    t1 = result_dict_n["t"][0]
    for p, t in zip(result_dict_n["p"], result_dict_n["t"]):
        s_list.append(t1 / t)
        e_list.append(t1 / (t * p))
    result_dict_n["s"] = s_list
    result_dict_n["e"] = e_list


def parse_results(k: int, net: int) -> dict:
    pattern = "results/result_{p}_{k}_{n}_net{net}.txt"
    # {net_size: {"p": [], "t": [], "edges": []}}
    results_dict = {}
    for n in [128, 256, 512]:
        results_dict[n] = {"p": [], "t": [], "edges": []}
        for p in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
            file_name = pattern.format(p=p, k=k, n=n, net=net)
            with open(file_name) as f:
                lines = [line.strip() for line in f.readlines()]
                if not lines:
                    print(f"{file_name} is empty")
            results_dict[n]["p"].append(p)
            results_dict[n]["t"].append(float(lines[0]))
            results_dict[n]["edges"].append(int(lines[1]))
        add_statistics(results_dict[n])
    return results_dict


def print_results(k: int, net: int, result_dict: dict) -> None:
    print(f"================ K = {k} ======= NET = {net} ================")
    print()
    print("\\begin{tabular}{ |c|c|c|c|c|c|}")
    print("\\hline")
    print("$p$   &   $size$   &   $T$   &   $S$   &   $E$   &   $edges$ \\\\")
    print("\\hline")
    for n in result_dict:
        print("\\hline")
        for i in range(len(result_dict[n]["p"])):
            print("{:3d} & {} & {:10.3f} & {:10.3f} & {:10.3f} & {:5d} \\\\".format(result_dict[n]["p"][i],
                                                                                    f"${n}\\times{n}$",
                                                                                    result_dict[n]["t"][i],
                                                                                    result_dict[n]["s"][i],
                                                                                    result_dict[n]["e"][i],
                                                                                    result_dict[n]["edges"][i]))
            print("\\hline")
    print("\\end{tabular}")
    print()
    print()


if __name__ == "__main__":
    visualize("out.txt")

    # for k in [16, 31, 64]:
    #     for net in [1, 2]:
    #         result = parse_results(k, net)
    #         print_results(k, net, result)
    #     make_graphics(result, k)
