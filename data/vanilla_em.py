#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.ticker import LogLocator


def update_Ex(Ex, Hz, Hy, cexe, cexh):
    Ex[:, 1:-1, 1:-1] = (cexe[:, 1:-1, 1:-1] *   Ex[:, 1:-1, 1:-1]
                      +  cexh[:, 1:-1, 1:-1] * ((Hz[:, 1:, 1:-1] - Hz[:, :-1, 1:-1])
                                             -  (Hy[:, 1:-1, 1:] - Hy[:, 1:-1, :-1])))

def update_Ey(Ey, Hx, Hz, ceye, ceyh):
    Ey[1:-1, :, 1:-1] = (ceye[1:-1, :, 1:-1] *   Ey[1:-1, :, 1:-1]
                      +  ceyh[1:-1, :, 1:-1] * ((Hx[1:-1, :, 1:] - Hx[1:-1, :, :-1])
                                             -  (Hz[1:, :, 1:-1] - Hz[:-1, :, 1:-1])))

def update_Ez(Ez, Hy, Hx, ceze, cezh):
    Ez[1:-1, 1:-1, :] = (ceze[1:-1, 1:-1, :] *   Ez[1:-1, 1:-1, :]
                      +  cezh[1:-1, 1:-1, :] * ((Hy[1:, 1:-1, :] - Hy[:-1, 1:-1, :])
                                             -  (Hx[1:-1, 1:, :] - Hx[1:-1, :-1, :])))

def update_Hx(Hx, Ey, Ez, chxh, chxe):
    Hx[:, :, :] = (chxh[:, :, :] * Hx[:, :, :]
                +  chxe[:, :, :] * ((Ey[:, :, 1:] - Ey[:, :, :-1])
                                 -  (Ez[:, 1:, :] - Ez[:, :-1, :])))

def update_Hy(Hy, Ez, Ex, chyh, chye):
    Hy[:, :, :] = (chyh[:, :, :] * Hy[:, :, :]
                +  chye[:, :, :] * ((Ez[1:, :, :] - Ez[:-1, :, :])
                                 -  (Ex[:, :, 1:] - Ex[:, :, :-1])))

def update_Hz(Hz, Ex, Ey, chzh, chze):
    Hz[:, :, :] = (chzh[:, :, :] * Hz[:, :, :]
                +  chze[:, :, :] * ((Ex[:, 1:, :] - Ex[:, :-1, :])
                                 -  (Ey[1:, :, :] - Ey[:-1, :, :])))


def ricker(n, cfl):
    loc = 0
    ppw = 25
    alpha = (np.pi * ((cfl * n - loc) / ppw - 1))**2
    return (1 - 2 * alpha) * np.exp(-alpha)


def plot(name, field, ax, fig):
    ax.clear()
    im = ax.contourf(field, levels=100)
    # fig.colorbar(im, ax=ax, format='{x:3.1e}', pad=0.01, shrink=0.5)
    # ax.set_aspect('equal')
    ax.set_title(f'{name}')

def main():
    nx = 101
    ny = 101
    nz = 101
    nt = 400
    cfl = 0.95 / np.sqrt(3)
    imp0 = 377.0

    Ex = np.zeros((nx - 1, ny, nz))
    Ey = np.zeros((nx, ny - 1, nz))
    Ez = np.zeros((nx, ny, nz - 1))
    Hx = np.zeros((nx, ny - 1, nz - 1))
    Hy = np.zeros((nx - 1, ny, nz - 1))
    Hz = np.zeros((nx - 1, ny - 1, nz))

    cexe = np.ones_like(Ex)
    ceye = np.ones_like(Ey)
    ceze = np.ones_like(Ez)
    cexh = np.full_like(Ex, cfl * imp0)
    ceyh = np.full_like(Ey, cfl * imp0)
    cezh = np.full_like(Ez, cfl * imp0)

    chxh = np.ones_like(Hx)
    chyh = np.ones_like(Hy)
    chzh = np.ones_like(Hz)
    chxe = np.full_like(Hx, cfl / imp0)
    chye = np.full_like(Hy, cfl / imp0)
    chze = np.full_like(Hz, cfl / imp0)

    save_step = 4
    count = 0

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8), layout='constrained', sharex=True, sharey=False)

    # fig, ax = plt.subplots()
    # ax.set_aspect('equal')
    # cax = make_axes_locatable(axes=ax).append_axes('right', '5%', '5%')

    for n in range(nt):
        print(f'Step {n}')
        update_Hx(Hx, Ey, Ez, chxh, chxe)
        update_Hy(Hy, Ez, Ex, chyh, chye)
        update_Hz(Hz, Ex, Ey, chzh, chze)

        Ez[nx // 2, ny // 2, nz // 2] += 10.0 * ricker(n, cfl)

        update_Ex(Ex, Hz, Hy, cexe, cexh)
        update_Ey(Ey, Hx, Hz, ceye, ceyh)
        update_Ez(Ez, Hy, Hx, ceze, cezh)


        if n % save_step == 0:
            plot('Ex', Ex[:, :, nz // 2], axes[0, 0], fig)
            plot('Ey', Ey[:, :, nz // 2], axes[0, 1], fig)
            plot('Ez', Ez[:, :, nz // 2], axes[0, 2], fig)
            plot('Hx', Hx[:, :, nz // 2], axes[1, 0], fig)
            plot('Hy', Hy[:, :, nz // 2], axes[1, 1], fig)
            plot('Hz', Hz[:, :, nz // 2], axes[1, 2], fig)
            fig.savefig(f'/home/cepheid/TriForce/game_engine/data/pngs/fields_{count:010d}.png')
            count += 1


if __name__ == '__main__':
    main()


