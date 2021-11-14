from math import log, ceil, sqrt
from typing import List


def compute_statistics(time_results: List[float], num_processes: List[int], n: int):
    print("==============={}^2={}===============".format(int(sqrt(n)), n))
    e_max_list = []
    s_p_list = []
    for p in num_processes:
        s_p = ceil(log(p, 2)) * (ceil(log(p, 2)) + 1) // 2
        e_max = 1. / (1. + log(p, n) * (log(p, 2) - 1.) / 2.)
        s_p_list.append(s_p)
        e_max_list.append(e_max)

    t_1 = time_results[0]
    print("  p          t             s      s_p       s_max          e          e_max       e_max percent")
    for p, t, s_p, e_max in zip(num_processes, time_results, s_p_list, e_max_list):
        s = t_1 / t
        e = s / p
        print("\\hline")
        print("{:3d} & {:10.3f} & {:10.3f} & {:3d} & {:10.3f} & {:10.3f} & {:10.3f} & {:10.3f} \\\\".format(p, t, s, s_p, e_max * p, e, e_max, e / e_max))

    print("\\hline")
    print()


if __name__ == "__main__":

    # {n: {num_processes: [time]}}
    result_dict = {}

    with open("results.txt", "r") as results_file:
        results = results_file.readlines()

    i = 0
    len_results = len(results)
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

    for n in sorted(result_dict):
        num_processes = sorted([n_p for n_p in result_dict[n]])
        time_results = [result_dict[n][n_p] for n_p in num_processes]
        compute_statistics(time_results, num_processes, n)

