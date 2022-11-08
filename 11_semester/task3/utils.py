from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os


def add_statistics(res_dict: dict) -> None:
    for key, val in res_dict.items():
        if key[0] == 1:
            res_dict[key]["s"] = 1
        else:
            res_dict[key]["s"] = res_dict[(1, key[1], key[2])]["t"] / res_dict[key]["t"]


def parse_results(results_dir: str) -> dict:
    results = {}

    for file_name in os.listdir(results_dir):
        if not file_name.startswith("out_"):
            continue
        _, num_p, N, L = file_name.split("_")
        L = "1" if L[:-4] == "1" else "\\pi"
        num_p = int(num_p)
        N = int(N)
        results[(num_p, N, L)] = {"t": [], "err": []}

        with open(os.path.join(results_dir, file_name), "r") as f:
            lines = [line.strip() for line in f.readlines()]
        for line in lines:
            _, _, error, time = line.split(" ")
            results[(num_p, N, L)]["t"].append(float(time))
            results[(num_p, N, L)]["err"].append(float(error))
        results[(num_p, N, L)]["t"] = np.mean(results[(num_p, N, L)]["t"])
        results[(num_p, N, L)]["err"] = np.mean(results[(num_p, N, L)]["err"])

    add_statistics(results)
    return results


def visualize(res_dict: dict, out_dir: str, L: str) -> None:
    text_dict = {"t": "Время выполнения", "s": "Ускорение"}
    font = {'family': 'normal',
            'weight': 'bold',
            'size': 24}

    matplotlib.rc('font', **font)
    for key in ("t", "s"):

        plt.figure(figsize=(15, 15))
        plt.title("Результаты ({}) для L={}".format(text_dict[key], L))
        plt.xlabel('Число процессов')
        plt.ylabel(text_dict[key])
        x = defaultdict(list)
        y = defaultdict(list)
        for key1 in sorted(res_dict.keys(), key=lambda x: (x[1], x[0])):
            x[key1[1]].append(key1[0])
            y[key1[1]].append(res_dict[key1][key])

        for key1 in y.keys():
            plt.plot(x[key1], y[key1], label=f"{key1} x {key1} x {key1}")
        loc = "upper right" if key == "t" else "upper left"
        plt.legend(loc=loc)
        plt.savefig(os.path.join(out_dir, f"{L}_{key}.jpg"))


if __name__ == "__main__":
    dir_out = "img"
    os.makedirs(dir_out, exist_ok=True)
    result = parse_results("polus.out")

    for L in ("1", "pi"):
        tex_l = "1" if L == "1" else "\\" + L
        print("\\begin{table}[H]")
        print("\\centering")
        print("\\begin{tabular}{|p{3cm}|p{2.3cm}|p{2.3cm}|p{2.7cm}|p{3.3cm}|}")
        print("\\hline")
        print("\\textbf{Число MPI-процессов} & \\textbf{Число точек сетки $N^3$} & \\textbf{Время решения $T$} & \\textbf{Ускорение $S$} & \\textbf{Погрешность $\\delta$ } \\\\")
        print("\\hline")
        filtered_result = {key: value for key, value in result.items() if key[2] == tex_l}
        for key in sorted(filtered_result.keys(), key=lambda x: (x[1], x[0])):
            print(f"{key[0]} & ${key[1]}^3$ & {filtered_result[key]['t']:10.3f} & {filtered_result[key]['s']:10.3f} & {filtered_result[key]['err']:.1E} \\\\")
        print("\\hline")
        print("\\end{tabular}")
        text = "{" + f"Результаты для Polus при $L={tex_l}$" + "}"
        print(f"\\caption{text}")
        print("\\label{table_" + L + "}")
        print("\\end{table}")
        print()
        print()
        visualize(filtered_result, dir_out, L)
