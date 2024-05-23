import os
import numpy as np
import matplotlib.pyplot as plt

base_path = os.path.join('..', 'output')
# 粒子配置字典
particle_configs = {
    'p': {
        'xlim': [0.0, 6.0],
        'plot_xlim': [0.0, 2.5],
        'ylim': [-1.6, 0.1],
        'file_prefix': 'proton_flow.dat',
        'file_exp_v1': 'v1_proton.dat',
        'file_exp_v2': 'v2_proton.dat',
        'label': r'p',
        'title': 'v1_AuAu_p_r.jpg',
        'folder_path': 'smash',
        'n_folders': 500,
        'output_file': 'v1_proton_sim.dat',
    },
    'd': {
        'xlim': [0.0, 6.0],
        'plot_xlim': [0.0, 1.7],
        'ylim': [-1.6, 0.1],
        'file_prefix': 'Deuteron_flow.dat',
        'file_exp_v1': 'v1_d.dat',
        'file_exp_v2': 'v2_d.dat',
        'label': r'd',
        'title': 'v1_AuAu_d_r.jpg',
        'folder_path': 'smash',
        'n_folders': 500,
        'output_file': 'v1_deuteron_sim.dat',
    },
    'he4': {
        'xlim': [0.0, 6.0],
        'plot_xlim': [0.0, 1.0],
        'ylim': [-1.6, 0.1],
        'file_prefix': 'Helium4_flow.dat',
        'file_exp_v1': 'v1_he4.dat',
        'file_exp_v2': 'v2_he4.dat',
        'label': r'$^4$He',
        'title': 'v1_AuAu_he4_r.jpg',
        'folder_path': 'smash',
        'n_folders': 160,
        'output_file': 'v1_he4_sim.dat',
    },
}

color_list = ['red', 'black', 'blue', 'green', 'purple', 'orange']
particle_colors = dict(zip(particle_configs.keys(), color_list))
legend_labels = {'exp': False, 'sim': False}
plt.figure(figsize=(9, 12))

# 遍历每种粒子
for particle, config in particle_configs.items():
    data_by_rapidity = {}
    experimental_data = {}
    n_folders = config['n_folders']
    output_path = os.path.join('..', 'output', config['output_file'])
    filename = config['file_prefix']
    experimental_v1_data_path = config['file_exp_v1']
    experimental_v2_data_path = config['file_exp_v2']
    plot_xlim = config['plot_xlim']

    with open(experimental_v1_data_path, 'r') as file:
        rapidity_ranges = ['-0.1<y<0', '-0.2<y<-0.1', '-0.3<y<-0.2',
                           '-0.4<y<-0.3']
        current_rapidity_index = 0

        for line in file:
            if line.strip() == '':
                current_rapidity_index += 1
            else:
                pt, v = map(float, line.split())
                if rapidity_ranges[
                    current_rapidity_index] not in experimental_data:
                    experimental_data[
                        rapidity_ranges[current_rapidity_index]] = []
                experimental_data[
                    rapidity_ranges[current_rapidity_index]].append(
                    (pt, v))

    # 读取数据
    for i in range(n_folders):
        folder_path = os.path.join(base_path, str(i))
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
                        pt, _, v1, _ = map(float, parts)
                        if pt not in data_by_rapidity[current_rapidity_range]:
                            data_by_rapidity[current_rapidity_range][pt] = []
                        data_by_rapidity[current_rapidity_range][pt].append(v1)

    # 计算统计数据
    particle_color = particle_colors[particle]
    stats_by_rapidity = {}
    for rapidity_range, pts_data in data_by_rapidity.items():
        stats_by_rapidity[rapidity_range] = {}
        for pt, v1s in pts_data.items():
            v1_mean = np.mean(v1s)
            v1_variance = np.var(v1s)
            stats_by_rapidity[rapidity_range][pt] = (v1_mean, v1_variance)

    # 绘图
    for rapidity_range, stats in stats_by_rapidity.items():
        if rapidity_range in experimental_data:
            multiplier = rapidity_multipliers[rapidity_range]
            scaling_factor = scaling_factors[particle]
            pts = np.array([pt / scaling_factor for pt in list(stats.keys())])
            means = np.array([mean * multiplier / scaling_factor for mean, _ in
                              stats.values()])
            variances = np.array([variance for _, variance in stats.values()])
            std_devs = np.sqrt(variances) / np.sqrt(n_folders) * multiplier
            pts_filtered, means_filtered = zip(
                *[(pt, mean) for pt, mean in zip(pts, means) if
                  plot_xlim[0] <= pt <= plot_xlim[1]])
            std_devs = np.array([std_dev for pt, std_dev in zip(pts, std_devs)
                                 if plot_xlim[0] <= pt <= plot_xlim[1]])
            if not legend_labels['sim']:
                plt.errorbar(pts_filtered, means_filtered, yerr=std_devs,
                             fmt='o',
                             label='Sim',
                             color=particle_color, capsize=5, markersize=8)
                means_upper = means_filtered + std_devs
                means_lower = means_filtered - std_devs
                plt.fill_between(pts_filtered, means_lower, means_upper,
                                 color=particle_color, alpha=0.5)
                legend_labels['sim'] = True
            else:
                plt.errorbar(pts_filtered, means_filtered, yerr=std_devs,
                             fmt='o',
                             color=particle_color, capsize=5, markersize=8)
                means_upper = means_filtered + std_devs
                means_lower = means_filtered - std_devs
                plt.fill_between(pts_filtered, means_lower, means_upper,
                                 color=particle_color, alpha=0.5)

    for rapidity_range, data in experimental_data.items():
        scaling_factor = scaling_factors[particle]
        data = np.array(data) / scaling_factor
        if not legend_labels['exp']:
            plt.plot(data[:, 0], data[:, 1],
                     label='Exp',
                     color=particle_color, marker='D', linestyle='--', markersize=8)
            legend_labels['exp'] = True
        else:
            plt.plot(data[:, 0], data[:, 1],
                     color=particle_color, marker='D', linestyle='--', markersize=8)

    with open(output_path, 'w') as file:
        for rapidity_range, pts_data in stats_by_rapidity.items():
            for pt, (mean, variance) in pts_data.items():
                file.write(f'{pt} {mean} {variance}\n')

for particle, color in particle_colors.items():
    plt.scatter([], [], color=color, label=particle_configs[particle]['label'],
                marker='o', s=100)

plt.title('AuAu 3GeV ($v_1$ vs $p_T$) 10%-40%', fontsize=28)
plt.xlabel('$p_T/A(GeV/c)$', fontsize=28)
plt.ylabel('$v_1/A$', fontsize=28)
plt.grid(True)
plt.legend(fontsize=16, loc='lower left')
plt.xlim(0.04, 2.8)
plt.ylim(-0.7, 0.1)
plt.text(2.0, -0.03, '-0.1<y<0', fontsize=16)
plt.text(2.0, -0.15, '-0.2<y<-0.1 (* 1.5)', fontsize=16)
plt.text(2.0, -0.3, '-0.3<y<-0.2 (* 2.0)', fontsize=16)
plt.text(2.0, -0.5, '-0.4<y<-0.3 (* 2.5)', fontsize=16)
plt.savefig('simp_smash_pdh2.jpg')
plt.show()
