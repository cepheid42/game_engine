import numpy as np
import matplotlib.pyplot as plt

nx, nz = 1500, 1500
xmin, xmax = -15.0e-6, 15.0e-6
zmin, zmax = -15.0e-6, 15.0e-6

# xs = np.linspace(xmin, xmax, nx, endpoint=True)
r = np.linspace(zmin, zmax, nz, endpoint=True)
# x, r = np.meshgrid(xs, rs)
x = -15.0e-6

w0 = 2.54796540e-6
wvl = 8.0e-7
E0 = 2.75e13

k = 2 * np.pi / wvl
xR = np.pi * w0**2 / wvl
wx = w0 * np.sqrt(1 + (x / xR)**2)
Rc = x * (1 + (xR / x)**2)
gouy = np.arctan(x / xR)
#
# print(k)
# print(x)
# print(xR)
# print(wx)
# print(Rc)
# print(gouy)

beam = (w0 / wx) * np.exp(-(r / wx)**2) * np.cos(k * x + 0.5 * k * r**2 / Rc - gouy)

coeffs = np.genfromtxt("/home/cepheid/TriForce/game_engine/data/coeffs.txt")

plt.plot(beam)
plt.plot(coeffs, '--')
plt.show()

# pulse = np.sin(np.pi * xs / 15.0e-6 + np.pi / 2)
# pulse[:375] = 0.0
# pulse[1124:] = 0.0
# beam = pulse * E0 * (w0 / wx) * np.exp(-r**2 / wx**2) * np.cos(k * x + 0.5 * k * r**2 / Rc - gouy)

# plt.contourf(x, r, beam, levels=100, cmap='jet')
# plt.show()

