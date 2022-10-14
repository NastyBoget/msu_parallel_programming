from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

esp2tex = {
    "1e-4": "$1.0 \\cdot 10^{-4}$",
    "2e-5": "$2.0 \\cdot 10^{-5}$",
    "8e-6": "$8.0 \\cdot 10^{-6}$",
    "3e-5": "$3.0 \\cdot 10^{-5}$",
    "5e-6": "$5.0 \\cdot 10^{-6}$",
    "1.5e-6": "$1.5 \\cdot 10^{-6}$"
}

tex2eps = {
    "$1.0 \\cdot 10^{-4}$": "1e-4",
    "$2.0 \\cdot 10^{-5}$": "2e-5",
    "$8.0 \\cdot 10^{-6}$": "8e-6",
    "$3.0 \\cdot 10^{-5}$": "3e-5",
    "$5.0 \\cdot 10^{-6}$": "5e-6",
    "$1.5 \\cdot 10^{-6}$": "1.5e-6"
}


def add_statistics(res_dict: dict) -> None:
    for key, val in res_dict.items():
        if key[0] == 1:
            res_dict[key]["s"] = 1
        else:
            res_dict[key]["s"] = res_dict[(1, key[1])]["t"] / res_dict[key]["t"]


def parse_results(results_dir: str) -> dict:
    results = {}

    for file_name in os.listdir(results_dir):
        if not file_name.startswith("out_"):
            continue
        _, num_p, eps = file_name.split("_")
        eps = esp2tex[eps[:-4]]
        num_p = int(num_p)
        results[(num_p, eps)] = {"t": [], "err": []}

        with open(os.path.join(results_dir, file_name), "r") as f:
            lines = [line.strip() for line in f.readlines()]
        for line in lines:
            int_val, err, points_num, time = line.split(" ")
            results[(num_p, eps)]["t"].append(float(time))
            results[(num_p, eps)]["err"].append(float(err))
        results[(num_p, eps)]["t"] = np.mean(results[(num_p, eps)]["t"])
        results[(num_p, eps)]["err"] = np.mean(results[(num_p, eps)]["err"])

    add_statistics(results)
    return results

def visualize(res_dict: dict, out_dir: str, prefix: str) -> None:
    text_dict = {"t": "Время выполнения", "s": "Ускорение"}
    font = {'family': 'normal',
            'weight': 'bold',
            'size': 24}

    matplotlib.rc('font', **font)
    for key in ("t", "s"):

        plt.figure(figsize=(15, 15))
        plt.title("Результаты ({}) для {}".format(text_dict[key], prefix))
        plt.xlabel('Число процессов')
        plt.ylabel(text_dict[key])
        x = defaultdict(list)
        y = defaultdict(list)
        for key1 in sorted(res_dict.keys(), key=lambda x: (x[1], x[0])):
            x[key1[1]].append(key1[0])
            y[key1[1]].append(res_dict[key1][key])

        for key1 in y.keys():
            plt.plot(x[key1], y[key1], label=tex2eps[key1])
        plt.legend(loc="upper right")
        plt.savefig(os.path.join(out_dir, f"{prefix}_{key}.jpg"))


if __name__ == "__main__":
    dir_out = "img"
    for cluster_name in ["polus", "bluegene"]:
        result = parse_results(f"{cluster_name}.out")

        print(f"======================{cluster_name}======================")
        print()
        print("\\begin{tabular}{|c|p{3cm}|p{3cm}|c|c|}")
        print("\\hline")
        print("\\textbf{Точность} $\\varepsilon $ & \\textbf{Число MPI-процессов} & \\textbf{Время работы программы} & \\textbf{Ускорение} & \\textbf{Ошибка} \\\\")
        print("\\hline")
        for key in sorted(result.keys(), key=lambda x: (x[1], x[0])):
            print(f"{key[1]} & {key[0]} & {result[key]['t']:10.3f} & {result[key]['s']:10.3f} & {result[key]['err']:.1E} \\\\")
        print("\\hline")
        print("\\end{tabular}")
        print()
        print()

        visualize(result, dir_out, cluster_name)
