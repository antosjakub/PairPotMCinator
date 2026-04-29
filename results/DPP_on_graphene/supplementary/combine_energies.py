


#####################################################################

def get_cycle_en_info(line):
    line = line.strip().replace(' |', ',')
    parts = line.split(',')
    assert len(parts) == 2
    en = parts[1].split(':')[1][1:]
    return float(en)

#line = 'cycle: 6 | en: -4.08049'
#print(get_cycle_en_info(line))

#####################################################################

def extract_energies(file):
    energies = []
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if '| en:' in line:
                en = get_cycle_en_info(line)
                energies.append(en)
    return energies

#####################################################################

def combine_energies(files):
    import numpy as np
    combined_energies = []
    n_files = len(files)
    for i in range(n_files):
        combined_energies.extend(extract_energies(files[i]))
    combined_energies = np.array(combined_energies)
    return combined_energies

#####################################################################

def write_energies(energies, file='combined_ens.txt'):
    with open(file, 'w') as f:
        n_energies = len(energies)
        for i in range(n_energies):
            f.write("\n")
            f.write(f"cycle: {i+1} | en: {energies[i]}")

#####################################################################

if __name__ == "__main__":

    filename_fn = lambda i: f'results_DTDPP_1x3_{i+1}/txt_files/cycle_min.txt'
    filenames = [filename_fn(i) for i in range(6)]

    energies = combine_energies(filenames)
    write_energies(energies, 'combined_ens.txt')