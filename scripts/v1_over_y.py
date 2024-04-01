import os
import numpy as np
import matplotlib.pyplot as plt

base_path = os.path.join('..', 'output', 'norandom')
folder_name = '10-40'
n_folders = 100
filename = 'proton_v1_3sub_plus.dat'
title = 'v1_plus_over_pt.jpg'


experimental_data_path = 'v1_over_y.dat'

experimental_data = []
with open(experimental_data_path, 'r') as file:
    for line in file:
        y, v2 = map(float, line.split())
        experimental_data.append((y, v2))

experimental_data = np.array(experimental_data)

data_by_rapidity = {}
for i in range(n_folders):
    folder_path = os.path.join(base_path, str(i))
    file_path = os.path.join(folder_path, filename)

    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()

        for line in lines:
            if line.startswith("Rapidity range"):
                parts = line.strip().split(',')
                rapidity_range_str = parts[0].split(':')[1].strip()
                y_values = [float(val) for val in rapidity_range_str.replace('<y<', ' ').split()]
                y_mid = np.mean(y_values)

                rapidity_range = f"{y_values[0]}<y<{y_values[1]}"
                if rapidity_range not in data_by_rapidity:
                    data_by_rapidity[rapidity_range] = {'x': [], 'y': []}
            elif line.startswith("flow over pt range"):
                v2_value = float(line.strip().split(',')[-1].split(':')[-1])
                data_by_rapidity[rapidity_range]['x'].append(y_mid)
                data_by_rapidity[rapidity_range]['y'].append(v2_value)


stats_by_rapidity = {}
for rapidity_range, xy_data in data_by_rapidity.items():
    x_mean = np.mean(xy_data['x'])
    y_mean = np.mean(xy_data['y'])
    y_variance = np.var(xy_data['y'])
    stats_by_rapidity[rapidity_range] = (x_mean, y_mean, y_variance)

plt.figure(figsize=(12, 8))

sim_color = 'blue'
exp_color = 'red'

for rapidity_range, (x_mean, y_mean, y_variance) in stats_by_rapidity.items():
    y_error = np.sqrt(y_variance)
    plt.errorbar(x_mean, y_mean, yerr=y_error, fmt='o', color=sim_color, label=f'Sim {rapidity_range}', capsize=5)

plt.plot(experimental_data[:, 0], experimental_data[:, 1], label='Exp Data', color=exp_color, marker='D', linestyle='--')

plt.title('Comparison of Simulation and Experimental Data ($v_1$ vs $y$)')
plt.xlabel('Rapidity (y)')
plt.ylabel('$v_1$')
plt.grid(True)
plt.legend()

plt.savefig(title)
plt.show()
