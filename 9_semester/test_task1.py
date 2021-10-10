import json
import os
import subprocess


if __name__ == "__main__":
    os.system("make task1")
    with open("test_data_task1.json") as f:
        test_data = json.load(f)
    for arr_size in test_data:
        output = subprocess.run(["./bsort", arr_size], capture_output=True, text=True)
        output = output.stdout.split('\n')
        n_comp, n_tact = output[-3], output[-2]
        assert(int(n_comp) <= test_data[arr_size][0])
        print(f"test passed for {arr_size} elements")
    os.remove("bsort")
