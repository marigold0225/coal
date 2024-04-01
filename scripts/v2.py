import os
import numpy as np
import matplotlib.pyplot as plt

# 设置基础路径和文件名
base_path = os.path.join('..', 'output', 'data')
folder_name = '10-40'
n_folders = 100
filename = 'proton_flow_10-40.dat'
title = 'v2_plus_10-40.jpg'
experimental_data_path = 'v2_proton.dat'

# 读取并解析实验数据文件
experimental_data = {}
with open(experimental_data_path, 'r') as file:
    # rapidity_ranges = ['-0.1<y<0', '-0.2<y<-0.1', '-0.3<y<-0.2', '-0.4<y<-0.3']
    rapidity_ranges = ['-0.1<y<0']
    current_rapidity_index = 0

    for line in file:
        if line.strip() == '':  # 新的快度区间块开始
            current_rapidity_index += 1
        else:
            pt, v2 = map(float, line.split())
            if rapidity_ranges[current_rapidity_index] not in experimental_data:
                experimental_data[rapidity_ranges[current_rapidity_index]] = []
            experimental_data[rapidity_ranges[current_rapidity_index]].append(
                (pt, v2))

data_by_rapidity = {}

for i in range(n_folders):
    folder_path = os.path.join(base_path, str(i),folder_name)
    file_path = os.path.join(folder_path, filename)

    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        if not lines:
            continue

        current_rapidity_range = None
        for line in lines:
            if line.startswith("Rapidity range"):
                range_str = line.strip().split(':')[1].split(',')[0].strip()
                current_rapidity_range = range_str
                if current_rapidity_range not in data_by_rapidity:
                    data_by_rapidity[current_rapidity_range] = {}
            elif current_rapidity_range:
                parts = line.split()
                if len(parts) == 4:
                    pt, _, _, v2 = map(float, parts)
                    if pt not in data_by_rapidity[current_rapidity_range]:
                        data_by_rapidity[current_rapidity_range][pt] = []
                    data_by_rapidity[current_rapidity_range][pt].append(v2)

stats_by_rapidity = {}
for rapidity_range, pts_data in data_by_rapidity.items():
    stats_by_rapidity[rapidity_range] = {}
    for pt, v2s in pts_data.items():
        v2_mean = np.mean(v2s)
        v2_variance = np.var(v2s)
        stats_by_rapidity[rapidity_range][pt] = (v2_mean, v2_variance)

# # 输出结果
# for rapidity_range, stats in stats_by_rapidity.items():
#     print(f"{rapidity_range}:")
#     for pt, (mean, variance) in stats.items():
#         print(f"  pt={pt}: Mean of v2 = {mean}, Variance of v2 = {variance}")

color_list = ['red', 'black', 'blue', 'green', 'purple']
colors = color_list[:len(experimental_data)]
color_map = dict(zip(experimental_data.keys(), colors))

plt.figure(figsize=(12, 8))

for rapidity_range, stats in stats_by_rapidity.items():
    if rapidity_range in experimental_data:  # 只绘制存在于实验数据中的快度区间
        pts = np.array(list(stats.keys()))
        means = np.array([mean for mean, _ in stats.values()])
        variances = np.array([variance for _, variance in stats.values()])
        std_devs = np.sqrt(variances)
        color = color_map[rapidity_range]
        plt.errorbar(pts, means, yerr=std_devs, fmt='o',
                     label=f'Sim {rapidity_range}', color=color, capsize=5,
                     markersize=8)
        # plt.plot(pts, means, marker='o', label=f'Sim {rapidity_range}',color=color)

for rapidity_range, data in experimental_data.items():
    data = np.array(data)
    color = color_map[rapidity_range]
    plt.plot(data[:, 0], data[:, 1], label=f'Exp {rapidity_range}',
             color='black', marker='D', linestyle='--')

plt.title('AuAu 3GeV ($v_1$ vs $p_T$)')
plt.xlabel('$p_T$')
plt.ylabel('$v_1$')
plt.xlim(0.4, 1.75)
# plt.ylim(-0.2, 0.2)
plt.grid(True)
plt.legend()
plt.savefig(title)
plt.show()
