"""
Imagine your 100-cycle simulation just finished and you realized that the energy minima are all listed somewhere in the all_min.txt file.
You would like to have a file that just list the minimum after each cycle, no other fuss.
Well this is where this script comes in.
Input: all_min.txt
Output: cycle_min_reconstructed.txt
"""


import sys

if len(sys.argv) > 1:
    file_in = sys.argv[1]
else:
    file_in = "results_DTDPP_1x3_1/txt_files/all_min.txt"

#####################################################################


def get_cycle_en_info(line):
    line = line.strip().replace(' |', ',')
    parts = line.split(',')
    assert len(parts) == 2
    en = parts[1].split(':')[1][1:]
    return float(en)

line = 'cycle: 1, cycle_min: -7.6331'
print(get_cycle_en_info(line))

#####################################################################

energies = []
with open(file_in, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if 'cycle_min' in line:
            en = get_cycle_en_info(line)
            energies.append(en)
print(energies)

#####################################################################

# each like should be smthg like this:
#cycle: 29 | en: -7.61326

import os

dirpath = os.path.dirname(file_in)
file_out = "cycle_min_reconstructed.txt"
file_out = os.path.join(dirpath, file_out)

n_cycles = len(energies)
with open(file_out, "w") as f:
    f.write('\n')
    for i in range(n_cycles):
        f.write(f'cycle: {i+1} | en: {energies[i]}\n')


