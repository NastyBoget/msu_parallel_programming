from math import log, ceil

time_results = [35.3468, 18.9226, 10.8192, 6.72132, 4.19626, 2.64018, 1.65177, 1.04559, 0.649666]
num_processes = [1, 2, 4, 8, 16, 32, 64, 128, 256]
n = 4096 * 4096


if __name__ == "__main__":
    e_max_list = []
    for p in num_processes:
        s_p = ceil(log(p, 2)) * (ceil(log(p, 2)) + 1) / 2.
        e_max = log(n, 2) / (log(n, 2) + s_p - log(p, 2))
        e_max_list.append(e_max)
    for p, e_max in zip(num_processes, e_max_list):
        print("p={} e_max={:.3f}".format(p, e_max))

    print()

    t_1 = time_results[0]
    e_list = []
    for p, t in zip(num_processes, time_results):
        s = t_1 / t
        e = s / p
        e_list.append(e)
        print("p={} s={:.3f} e={:.3f}".format(p, s, e))

    print()

    for p, e, e_max in zip(num_processes, e_list, e_max_list):
        if e < 0.5 * e_max:
            print("Too low effectiveness for {} processes".format(p))
