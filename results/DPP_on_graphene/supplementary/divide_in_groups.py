import numpy as np

"""
labels:

0: random
1: row tile (c)
2: column-row (b)
3: fork (a)
4: glob min (d)

run 1:
0303133213
run 2:
2303133013
run 3:
2303032013
run 4:
3303032013
run 5:
2333132043
run 6:
23331322033233331321102132322210321031232223211113
"""


folder = '/home/antos_j/molecular_configuration_UK/results/DPP_on_graphene/supplementary/'
filename = 'combined_trajs.xyz'
filename = f'{folder}/{filename}'
labels = [
'0303133213',
'2303133013',
'2303032013',
'3303032013',
'2333132043',
'23331322433233331321102132322214321431232223211113'
]
labels = ''.join(labels)
labels = list(map(lambda x: int(x), labels))

# input
# - filename
# - labels

assert min(labels) == 0


from combine_trajectories import read_trajs

trajs = read_trajs(filename)
n_trajs = len(trajs)

print(f"Number of trajs = {n_trajs}")

traj_groups = []
traj_groups_metadata = []
labels_unique = np.unique(labels)
for l in labels_unique:
    traj_groups.append([])
    traj_groups_metadata.append([])

for i, (traj, label) in enumerate(zip(trajs, labels)):
    #traj_groups[label].append('\n'.join(traj))
    traj_groups[label].extend(traj)
    traj_groups_metadata[label].append(traj[1])


folder = '3_groups'
for k in labels_unique:
    with open(f'{folder}/superstructs_label={k}.xyz', "w") as f:
        f.writelines('\n'.join(traj_groups[k]))
with open(f'{folder}/superstructs_groups_metadata.txt', "w") as f:
    for k in labels_unique:
        group_metadata = traj_groups_metadata[k]
        n_trajs_in_group = len(group_metadata)
        f.write(f'-------- LABEL = {k} --------\n')
        for i in range(n_trajs_in_group):
            f.write(f'{i+1}: ({group_metadata[i]})\n')
        f.write(f'---------------------------\n')
