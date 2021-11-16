from math import log, ceil, sqrt
from typing import List
import matplotlib.pyplot as plt
import numpy as np


def compute_statistics(time_results: List[float], num_processes: List[int], n: int) -> dict:
    print("==============={}^2={}===============".format(int(sqrt(n)), n))
    e_max_list = []
    s_p_list = []
    for p in num_processes:
        s_p = ceil(log(p, 2)) * (ceil(log(p, 2)) + 1) // 2
        e_max = 1. / (1. + log(p, n) * (log(p, 2) - 1.) / 2.)
        s_p_list.append(s_p)
        e_max_list.append(e_max)

    t_1 = time_results[0]
    print("T1/[N*log(N)] = {}".format(t_1/(n*log(n, 2))))
    result = {"p": [], "t": [], "s": [], "s_p": [], "s_max": [], "e": [], "e_max": []}
    print("  p          t             s      s_p       s_max          e          e_max       e_max percent")
    for p, t, s_p, e_max in zip(num_processes, time_results, s_p_list, e_max_list):
        s = t_1 / t
        e = s / p
        print("\\hline")
        print("{:3d} & {:10.3f} & {:10.3f} & {:3d} & {:10.3f}"
              " & {:10.3f} & {:10.3f} & {:10.3f} \\\\".format(p, t, s, s_p, e_max * p, e, e_max, e / e_max))
        result["p"].append(p)
        result["t"].append(t)
        result["s"].append(s)
        result["s_p"].append(s_p)
        result["s_max"].append(e_max * p)
        result["e"].append(e)
        result["e_max"].append(e_max)

    print("\\hline")
    print()
    return result


def parse_results(results_path: str) -> dict:
    # {n: {num_processes: [time]}}
    result_dict = {}

    with open(results_path, "r") as results_file:
        results = results_file.readlines()

    i = 0
    while i < len(results) - 1:
        info = results[i].strip().split(" ")
        if len(info) != 3:
            i += 1
            continue
        n = int(info[0]) * int(info[1])
        n_processes = int(info[2])
        info = results[i + 1].strip().split(" ")
        if len(info) != 1:
            i += 1
            continue
        time_result = float(info[0])
        if n in result_dict:
            result_dict[n][n_processes] = time_result
        else:
            result_dict[n] = {n_processes: time_result}
        i += 2
    return result_dict


def make_graphics(statistics: dict) -> None:
    text_dict = {"t": "Время выполнения", "s": "Ускорение", "e": "Эффективность"}
    color_dict = {"t": '-bD', "s": '-gD', "e": '-rD'}

    for key in ("t", "s", "e"):
        fig, axes = plt.subplots(2, 2)
        fig.set_figheight(15)
        fig.set_figwidth(15)
        for k, n in enumerate(statistics):
            i, j = k // 2, k % 2
            x = statistics[n]["p"]
            y = np.array(statistics[n][key])

            axes[i, j].set(title="{}: результаты для сетки {}x{}".format(text_dict[key], int(sqrt(n)), int(sqrt(n))),
                           xlabel='Число процессов', ylabel=text_dict[key])

            axes[i, j].plot(x, y, color_dict[key], markevery=np.arange(len(y)))
            for i_x in range(len(x)):
                axes[i, j].annotate(str((x[i_x], round(y[i_x], 3))), (x[i_x], y[i_x]))

        fig.savefig('{}.jpg'.format(key))


if __name__ == "__main__":
    result_dict = parse_results("results.txt")
    statistics = {}
    for n in sorted(result_dict):
        num_processes = sorted([n_p for n_p in result_dict[n]])
        time_results = [result_dict[n][n_p] for n_p in num_processes]
        statistics[n] = compute_statistics(time_results, num_processes, n)

