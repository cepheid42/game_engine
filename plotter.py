import numpy as np
import matplotlib.pyplot as plt


data_path = '/home/cepheid/TriForce/game_engine/data'

nx = 100
ny = 200
nz = 100
nt = 200

# times = 1e-9 * np.array([1673642100, 1629794566, 3306171535, 1750413224, 1610065498, 3363739950])
#
# ave_times = times / nt
#
# print(f'Fancy E update: total = {times[0]}, per step = {ave_times[0]}')
# print(f'Fancy B update: total = {times[1]}, per step = {ave_times[1]}')
# print(f'Fancy Total: total = {times[2]}, per step = {ave_times[2]}')
# print(f'Regular E update: total = {times[3]}, per step = {ave_times[3]}')
# print(f'Regular B update: total = {times[4]}, per step = {ave_times[4]}')
# print(f'Regular Total: total = {times[5]}, per step = {ave_times[5]}')

save_step = 2
num_files = nt // save_step

for i in range(num_files):
    filename = f'/ez_{i:08d}.csv'
    print(f'Processing file: {filename[1:]}')
    data = np.loadtxt(data_path + filename, dtype=np.float64, delimiter=',').reshape((nx, ny))
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.contourf(data, levels=50)
    plt.colorbar(im)
    plt.savefig(data_path + f'/pngs' + filename[:-3] + 'png')

    plt.clf()
    plt.close(fig)
