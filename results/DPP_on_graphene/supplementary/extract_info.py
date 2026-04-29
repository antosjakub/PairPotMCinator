
import sys

if len(sys.argv) > 1:
    file = sys.argv[1]
    pass
else:
    file = "results_DTDPP_1x3_3/txt_files/all_min.txt"


n_cycles = 10
n_loops_stage_1 = 200

#####################################################################

line = 'cycle: 1, stage: 1, loop: 35 | en: -4.08049'

def get_cycle_stage_loop_en_info(line):
    line = line.strip().replace(' |', ',')
    parts = line.split(',')
    assert len(parts) == 4
    info = []
    for i in range(4):
        number = parts[i].split(':')[1][1:]
        if i != 3:
            info.append(int(number))
        else:
            info.append(float(number))
    return info

print(get_cycle_stage_loop_en_info(line))

#####################################################################

def extract_all(filename):
    info = []
    ens = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'cycle' in line and 'stage' in line and 'loop' in line:
                all = get_cycle_stage_loop_en_info(line)
                info.append(all[:-1])
                ens.append(all[-1])

    return info, ens

info, ens = extract_all(file)
print(len(info))
print(len(ens))
print(info[:4])
print(ens[:4])

print(20*'#')

import numpy as np
info_np = np.array(info)
ens_np = np.array(ens)

info_cycle = info_np[:,0]
info_stage = info_np[:,1]
info_loop = info_np[:,2]

filter = (info_cycle == 3) & (info_stage == 1)

info = info_np[filter]
ens = ens_np[filter]
for i in range(len(ens)):
    print(info[i], "en:", ens[i])




### extract en min for each cycle
print(10*'===---')
l = ''
for i in range(10):
    filter = (info_cycle == i+1)
    ens = ens_np[filter]
    en_min = ens[-1]
    print(i+1, en_min)
    l += str(round(en_min,3)) + ' & '
print(l)
print(10*'===---')



##############################################################################################################

from matplotlib import pyplot as plt

#for i in range(10):
#    filter = (info_cycle == i+1) & (info_stage == 1)
#    ens = ens_np[filter]
#    plt.plot(ens)
#plt.savefig('ens.png')

#x_loops = info[:,2]
#x_loops_log = np.log(x_loops)
#
#plt.scatter(x_loops_log, ens)
#plt.xscale('log')
#plt.xlim(1, 200)
#
#plt.xlabel("X (log scale)")
#plt.ylabel("Y")
#plt.savefig('ens.png')

#x = np.linspace(1, 200, 200)
#y = np.linspace(-10, 1000, 200)

for i in range(10):
    filter = (info_cycle == i+1) & (info_stage == 1)
    ens = ens_np[filter]
    info = info_np[filter]
    x = info[:,2]
    y = ens

##############################################################################################################
############### ENERGY VS LOOP NUMBER LOG LOG PLOT ###########################################################
##############################################################################################################

import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# Example data (REPLACE with your real data)
# x_data and y_data must be lists of arrays
# --------------------------------------------------
np.random.seed(0)

x_data = []
y_data = []

for i in range(n_cycles):
    #x = np.linspace(0, 600, 300)
    #y = 50*np.log1p(x) - 20 + i*15 + np.random.normal(0, 5, x.size)
    filter = (info_cycle == i+1)
    y = ens_np[filter]
    print(y.max())
    x1 = info_np[filter & (info_stage == 1)][:,2]
    x2 = info_np[filter & (info_stage == 2)][:,2] + n_loops_stage_1
    x = np.concatenate((x1,x2))
    x_data.append(x)
    y_data.append(y)

# --------------------------------------------------
# Plot styling (scientific & minimal)
# --------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.fontsize": 9
})

fig, ax = plt.subplots(figsize=(7.2, 4.5))

# Colorblind-safe default cycle
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['-', '--', '-.', ':']

# --------------------------------------------------
# Plot each curve
# --------------------------------------------------
for i, (x, y) in enumerate(zip(x_data, y_data)):
    # Log x-axis cannot include x <= 0
    mask = x > 0
    ax.plot(
        x[mask],
        y[mask],
        linewidth=1.4,
        linestyle=linestyles[i % len(linestyles)],
        color=colors[i % len(colors)],
        label=f"cycle {i+1}"
    )

# --------------------------------------------------
# Axes scaling and limits
# --------------------------------------------------
ax.set_xscale('log')
ax.set_yscale('symlog', linthresh=1)

ax.set_xlim(1, 600)
ax.set_ylim(-10, 1500)

ax.set_xlabel("loop")
ax.set_ylabel(r"$E$ [eV]")

# --------------------------------------------------
# Grid and legend
# --------------------------------------------------
ax.grid(True, which='both', linestyle=':', linewidth=0.6, alpha=0.7)

ax.legend(
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    frameon=False
)

plt.tight_layout()
plt.savefig('ens.png')