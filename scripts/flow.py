import os
import numpy as np
import matplotlib.pyplot as plt

# Global constants
BASE_PATH = os.path.join('..', 'output', 'smash')
RAPIDITY_MULTIPLIERS = {
    '-0.1<y<0': 1, '-0.2<y<-0.1': 1.5, '-0.3<y<-0.2': 2.0, '-0.4<y<-0.3': 2.5
}
SCALING_FACTORS = {
    'p': 1, 'd': 2, 'he4': 4, 'be8': 8
}
legend_labels = {'exp': False, 'sim': False}
particle_configs = {
    'p': {
        'xlim_pt': [0.0, 6.0],
        'xlim_v1': [0.0, 2.5],
        'xlim_v2': [0.0, 2.0],
        'ylim_v1': [-1.6, 0.1],
        'ylim_v2': [-1.6, 0.1],
        'file_prefix': 'proton_flow.dat',
        'file_exp_v1': 'v1_proton.dat',
        'file_exp_v2': 'v2_proton.dat',
        'file_exp_v1_over_y': 'v1_over_y.dat',
        'file_exp_v2_over_y': 'v2_over_y.dat',
        'label': 'p',
        'title': 'v1_AuAu_p_r.jpg',
        'n_folders': 500,
        'output_file': 'v1_proton_sim.dat',
    },
    'd': {
        'xlim_pt': [0.0, 6.0],
        'xlim_v1': [0.0, 1.7],
        'xlim_v2': [0.0, 2.0],
        'ylim_v1': [-1.6, 0.1],
        'ylim_v2': [-1.6, 0.1],
        'file_prefix': 'Deuteron_flow.dat',
        'file_exp_v1': 'v1_d.dat',
        'file_exp_v2': 'v2_d.dat',
        'file_exp_v1_over_y': 'v1_over_y.dat',
        'file_exp_v2_over_y': 'v2_over_y.dat',
        'label': 'd',
        'title': 'v1_AuAu_d_r.jpg',
        'n_folders': 500,
        'output_file': 'v1_deuteron_sim.dat',
    },
    'he4': {
        'xlim_pt': [0.0, 6.0],
        'xlim_v1': [0.0, 1.0],
        'xlim_v2': [0.0, 2.0],
        'ylim_v1': [-1.6, 0.1],
        'ylim_v2': [-1.6, 0.1],
        'file_prefix': 'Helium4_flow.dat',
        'file_exp_v1': 'v1_he4.dat',
        'file_exp_v2': 'v2_he4.dat',
        'file_exp_v1_over_y': 'v1_over_y.dat',
        'file_exp_v2_over_y': 'v2_over_y.dat',
        'label': '$^4$He',
        'title': 'v1_AuAu_he4_r.jpg',
        'n_folders': 160,
        'output_file': 'v1_he4_sim.dat',
    },
    'be8': {
        'xlim_pt': [0.0, 6.0],
        'xlim_v1': [0.0, 0.5],
        'xlim_v2': [0.0, 2.0],
        'ylim_v1': [-1.6, 0.1],
        'ylim_v2': [-1.6, 0.1],
        'file_prefix': 'Be8_flow.dat',
        # 'file_exp_v1': 'v1_be8.dat',
        # 'file_exp_v2': 'v2_be8.dat',
        # 'file_exp_v1_over_y': 'v1_over_y.dat',
        # 'file_exp_v2_over_y': 'v2_over_y.dat',
        'label': '$^8$Be',
        'title': 'v1_AuAu_be8_r.jpg',
        'n_folders': 1,
        'output_file': 'v1_be8_sim.dat',
    },
}
color_list = ['red', 'black', 'blue', 'green', 'purple', 'orange']
particle_colors = dict(zip(particle_configs.keys(), color_list))


def load_exp_flow_data(filepath):
    exp_data = {}
    with open(filepath, 'r') as file:
        rapidity_ranges = ['-0.1<y<0', '-0.2<y<-0.1', '-0.3<y<-0.2',
                           '-0.4<y<-0.3']
        current_rapidity_index = 0

        for line in file:
            if line.strip() == '':
                current_rapidity_index += 1
            else:
                pt, v = map(float, line.split())
                if rapidity_ranges[current_rapidity_index] not in exp_data:
                    exp_data[
                        rapidity_ranges[current_rapidity_index]] = []
                exp_data[
                    rapidity_ranges[current_rapidity_index]].append((pt, v))
    return exp_data


def load_exp_flow_over_y(particle, filepath):
    exp_data = {particle: []}
    particle_section = False

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line == f'{particle}:':
                particle_section = True
                continue
            if line.endswith(':') and particle_section:
                break
            if particle_section:
                y, v = map(float, line.split())
                exp_data[particle].append((y, v))

    return exp_data


def load_simulation_data(base_path, n_folders, filename, column_index):
    column_map = {
        'yield': 1, 'v1': 2, 'v2': 3, 'v1_over_y': 'v1', 'v2_over_y': 'v2'
    }
    data_by_rapidity = {}
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
                if column_index in ['yield', 'v1', 'v2']:
                    if line.startswith("Rapidity range"):
                        range_str = line.strip().split(':')[1].split(',')[
                            0].strip()
                        current_rapidity_range = range_str
                        if current_rapidity_range not in data_by_rapidity:
                            data_by_rapidity[current_rapidity_range] = {}
                    elif current_rapidity_range:
                        parts = line.split()
                        if len(parts) == 4:
                            pt = float(parts[0])
                            v = float(parts[column_map[column_index]])
                            if pt not in data_by_rapidity[
                                current_rapidity_range]:
                                data_by_rapidity[current_rapidity_range][
                                    pt] = []
                            data_by_rapidity[current_rapidity_range][pt].append(
                                v)

                elif column_index in ['v1_over_y', 'v2_over_y']:
                    if line.startswith("Rapidity range"):
                        range_str = line.strip().split(':')[1].split(',')[
                            0].strip()
                        y_values = [float(y) for y in
                                    range_str.replace('<y<', ' ').split()]
                        y_mid = np.mean(y_values)
                        current_rapidity_range = range_str
                        if current_rapidity_range not in data_by_rapidity:
                            data_by_rapidity[current_rapidity_range] = {'x': [],
                                                                        'y': []}
                    elif current_rapidity_range:
                        if line.startswith("flow over pt range"):
                            v_values_str = line.strip().split(',')
                            v_value = float([part.split(':')[1] for part in
                                             v_values_str if
                                             column_map[column_index] in part][
                                                0].strip())
                            data_by_rapidity[current_rapidity_range][
                                'x'].append(y_mid)
                            data_by_rapidity[current_rapidity_range][
                                'y'].append(v_value)
    return data_by_rapidity


def process_data(data_by_rapidity):
    stats_by_rapidity = {}
    for rapidity_range, data in data_by_rapidity.items():
        if 'x' in data:
            x_values = data['x']
            y_values = data['y']
            y_mean = np.mean(y_values)
            y_variance = np.var(y_values)
            stats_by_rapidity[rapidity_range] = (x_values, y_mean, y_variance)
        else:
            stats_by_rapidity[rapidity_range] = {}
            for pt, v in data.items():
                v_mean = np.mean(v)
                v_variance = np.var(v)
                stats_by_rapidity[rapidity_range][pt] = (v_mean, v_variance)
    return stats_by_rapidity


def plot_flow(stats_by_rapidity, experimental_data, n_folders, plot_xlim,
              particle_color, particle_i, label):
    for rapidity_range, stats in stats_by_rapidity.items():
        # if rapidity_range in experimental_data:
        multiplier = RAPIDITY_MULTIPLIERS[rapidity_range]
        pts = np.array(
            [pt / SCALING_FACTORS[particle_i] for pt in list(stats.keys())])
        means = np.array(
            [mean * multiplier / SCALING_FACTORS[particle_i] for mean, _ in
             stats.values()])
        variance = np.array([variance for _, variance in stats.values()])
        std_devs = np.sqrt(variance) / np.sqrt(n_folders) * multiplier
        pts_filtered, means_filtered = zip(*[(pt, mean) for
                                             pt, mean in zip(pts, means)
                                             if plot_xlim[0] <= pt <=
                                             plot_xlim[1]])
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

    if experimental_data:
        for rapidity_range, data in experimental_data.items():
            scaling_factor = SCALING_FACTORS[particle_i]
            data = np.array(data) / scaling_factor
            if not legend_labels['exp']:
                plt.plot(data[:, 0], data[:, 1],
                         label='Exp',
                         color=particle_color, marker='D', linestyle='--',
                         markersize=8)
                legend_labels['exp'] = True
            else:
                plt.plot(data[:, 0], data[:, 1],
                         color=particle_color, marker='D', linestyle='--',
                         markersize=8)


def plot_flow_over_y(particle, stats_by_rapidity, experimental_data, n_folders,
                     particle_color, label):
    for rapidity_range, stats in stats_by_rapidity.items():
        x_values = np.mean(stats['x'])
        y_values = np.mean(stats['y'])
        y_variance = np.var(stats['y'])
        y_error = np.sqrt(y_variance) / np.sqrt(n_folders)
        if not legend_labels['sim']:
            plt.errorbar(x_values, y_values, yerr=y_error, fmt='o',
                         color=particle_color,
                         label='Sim', capsize=5, markersize=8)
            legend_labels['sim'] = True
        else:
            plt.errorbar(x_values, y_values, yerr=y_error, fmt='o',
                         color=particle_color, capsize=5, markersize=8)

    if experimental_data:
        if particle in experimental_data:
            x, y = zip(*experimental_data[particle])
            if not legend_labels['exp']:
                plt.plot(x, y, label='Exp Data', color=particle_color,
                         marker='D',
                         linestyle='--')
                legend_labels['exp'] = True
            else:
                plt.plot(x, y, color=particle_color, marker='D', linestyle='--')


def flow_over_y(particles, plot_type):
    plt.figure(figsize=(12, 8))
    plt.rcParams['font.size'] = 18
    plt.rcParams['font.family'] = 'Times New Roman'

    for particle_i in particles:
        config = particle_configs.get(particle_i, {})
        exp_data_key = f'file_exp_{plot_type}'
        exp_data = None
        if exp_data_key in config:
            try:
                exp_data = load_exp_flow_over_y(particle_i,
                                                config[exp_data_key])
            except FileNotFoundError:
                print(f"File not found: {config[f'{plot_type}.dat']}")
        simulation_data = load_simulation_data(BASE_PATH, config['n_folders'],
                                               config['file_prefix'],
                                               plot_type)
        plot_flow_over_y(particle_i, simulation_data, exp_data,
                         config['n_folders'],
                         particle_colors[particle_i], config['label'])

    for particle, color in particle_colors.items():
        plt.scatter([], [], color=color,
                    label=particle_configs[particle]['label'],
                    marker='o', s=100)
    plt.title('AuAu 3GeV ($v_1$ vs $y$) 10%-40%', fontsize=28)
    plt.xlabel('Rapidity (y)', fontsize=28)
    plt.ylabel('$v_1$', fontsize=28)
    plt.xlim(-0.5, 0.0)
    # plt.ylim(-0.6, 0.05)
    plt.ylim(-0.05, 0.2)
    plt.grid(True)
    plt.legend(fontsize=16)
    plt.savefig(f'{plot_type}.jpg')
    plt.show()


def flow(particles, plot_type):
    plt.figure(figsize=(9, 12))
    plt.rcParams['font.size'] = 18
    plt.rcParams['font.family'] = 'Times New Roman'

    for particle_i in particles:
        config = particle_configs.get(particle_i, {})
        exp_data_key = f'file_exp_{plot_type}'
        exp_data = None
        if exp_data_key in config:
            exp_data_path = config[exp_data_key]
            try:
                exp_data = load_exp_flow_data(exp_data_path)
            except FileNotFoundError:
                print(f"File not found: {exp_data_path}")
        simulation_data = load_simulation_data(BASE_PATH, config['n_folders'],
                                               config['file_prefix'],
                                               plot_type)
        stats_by_rapidity = process_data(simulation_data)
        plot_flow(stats_by_rapidity, exp_data, config['n_folders'],
                  config[f'xlim_{plot_type}'],
                  particle_colors[particle_i], particle_i, config['label'])

    for particle, color in particle_colors.items():
        plt.scatter([], [], color=color,
                    label=particle_configs[particle]['label'],
                    marker='o', s=100)
    plt.title('AuAu 3GeV ($v_1$ vs $p_T$) 10%-40%', fontsize=28)
    plt.xlabel('$p_T/A(GeV/c)$', fontsize=28)
    plt.ylabel('$v_1/A$', fontsize=28)
    plt.grid(True)
    plt.legend(fontsize=16, loc='lower left')
    # plt.xlim(0.0, 2.0)
    # plt.ylim(-0.1,0.08)
    plt.xlim(0.04, 2.8)
    plt.ylim(-0.7, 0.1)
    plt.text(2.0, -0.03, '-0.1<y<0', fontsize=18)
    plt.text(2.0, -0.15, '-0.2<y<-0.1 (* 1.5)', fontsize=18)
    plt.text(2.0, -0.3, '-0.3<y<-0.2 (* 2.0)', fontsize=18)
    plt.text(2.0, -0.5, '-0.4<y<-0.3 (* 2.5)', fontsize=18)
    plt.savefig(f'{plot_type}.jpg')
    plt.show()


flow(['p', 'd', 'he4', 'be8'], 'v1')
# flow_over_y(['p', 'd', 'he4'], 'v2_over_y')
