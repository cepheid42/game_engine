import numpy as np

import matplotlib.pyplot as plt

nt = 400
step = 4

data_path = '/home/cepheid/TriForce/game_engine/data'

for n in range(0, nt // step):
    print(f'Step {n}')
    # file1 = data_path + f'_slim/Ez_{n:04d}.csv'
    file2 = data_path + f'/Ez_{n:04d}.csv'

    # data1 = np.genfromtxt(file1, dtype=np.float64, delimiter=',')
    data2 = np.genfromtxt(file2, dtype=np.float64, delimiter=',')

    fig, ax = plt.subplots()
    # ax.contourf(data)
    # ax.plot(data1, label='data1')
    ax.plot(data2, label='data2')
    ax.set_ylim([-1, 1])
    ax.legend()
    plt.savefig(data_path + f'/pngs/ez_{n:04d}.png')
    plt.clf()
    plt.close(fig)
