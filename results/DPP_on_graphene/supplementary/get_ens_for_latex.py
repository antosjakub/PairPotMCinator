
import sys
filename = sys.argv[1]


import numpy as np
from combine_energies import extract_energies
energies = extract_energies(filename)

string = ""
for i in range(len(energies)):
    string += f"{energies[i]:.3f} & "

print(string)