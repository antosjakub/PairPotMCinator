"""
Combine multiple .xyz files into a single one.
"""



def read_trajs(filename):
    """
    read each trajectory as a list of lines
    """
    trajs = []
    with open(filename, "r") as f:
        while True:
            traj_single = []
            try:
                # number of atoms
                line = f.readline().strip()
                traj_single.append(line)
                n_atoms = int(line)
                # structure name / comment
                line = f.readline().strip()
                traj_single.append(line)
            except:
                break

            for _ in range(n_atoms):
                line = f.readline().strip()
                traj_single.append(line)
                if not line:
                    # skip empty lines
                    continue

            trajs.append(traj_single)
    return trajs

def combine_trajs_into_one_file(filenames, prepend_traj_name_fn=lambda fi: f"traj {fi+1}, "):
    combined_trajs = []
    n_files = len(filenames)
    for i in range(n_files):
        trajs = read_trajs(filenames[i])
        n_trajs = len(trajs)
        for j in range(n_trajs):
            trajs[j][1] = prepend_traj_name_fn(i) + " " + trajs[j][1]
            combined_trajs.extend(['\n'.join(trajs[j])])
    return combined_trajs

def write_trajs(trajs, filename='combined_trajs.xyz'):
    with open(filename, "w") as f:
        f.writelines('\n'.join(trajs))


if __name__ == "__main__":

    # get the files
    folder = '/home/antos_j/molecular_configuration_UK/results/DPP_on_graphene/supplementary/'
    filename_fn = lambda i: f'{folder}/results_DTDPP_1x3_{i}/xyz_files/cycle_min_periodic.xyz'
    filenames = [filename_fn(i+1) for i in range(6)]

    prepend_traj_name_fn = lambda fi: f'run: {fi+1},' # fi = file index
    combined_trajs = combine_trajs_into_one_file(filenames, prepend_traj_name_fn)
    print(f"Total number of trajectories = {len(combined_trajs)}")
    write_trajs(combined_trajs, f"{folder}/combined_trajs.xyz")