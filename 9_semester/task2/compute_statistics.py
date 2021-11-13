from math import log, ceil
from typing import List

# time_results = [152.056, 81.7605, 47.428, 28.7055, 17.748, 11.0372, 6.84153, 4.3077, 2.6856]
# num_processes = [1, 2, 4, 8, 16, 32, 64, 128, 256]
# n = 8192 * 8192


def compute_statistics(time_results: List[float], num_processes: List[int], n: int):
    print("==============={}===============".format(n, n))
    e_max_list = []
    s_p_list = []
    for p in num_processes:
        s_p = ceil(log(p, 2)) * (ceil(log(p, 2)) + 1) // 2
        e_max = 1. / (1. + log(p, n) * (log(p, 2) - 1.) / 2.)
        s_p_list.append(s_p)
        e_max_list.append(e_max)
    for p, s_p, e_max in zip(num_processes, s_p_list, e_max_list):
        print("p={} s_p={} s_max={:.3f} e_max={:.3f}".format(p, s_p, e_max * p, e_max))

    print()

    t_1 = time_results[0]
    for p, t, e_max in zip(num_processes, time_results, e_max_list):
        s = t_1 / t
        e = s / p
        print("p={} t={:.3f} s={:.3f} e={:.3f} e_max percent={:.3f}".format(p, t, s, e, e / e_max))

    print()


if __name__ == "__main__":

    # {n: {num_processes: [time]}}
    result_dict = {}

    with open("results1.txt", "r") as results_file:
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

