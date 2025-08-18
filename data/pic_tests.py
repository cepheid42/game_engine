import numpy as np

xmin, xmax = 0, 10
nx = 11
dx = (xmax - xmin) / (nx - 1)

xold = np.array([(i * dx + 0.25) for i in range(10)])
xnew = np.copy(xold)
xnew[0] += 0.5
xnew[1] += 1
xnew[2] += 0.75
xnew[3] -= 0.5
xnew[4] -= 1
xnew[5] -= 0.75
xnew[6] += 0.5
xnew[7] += 0.9
xnew[8] += 0.1
xnew[9] -= 0.33

xold /= dx
xnew /= dx

print(xold)
print(xnew)

# xeven1 = np.floor(xold)
# xeven2 = np.floor(xold + 1)
#
# xodd1 = np.floor(xold + 0.5) - 0.5
# xodd2 = np.floor(xold + 0.5) + 0.5

# print(xeven1)
# print(xeven2)

for xo, xn in zip(xold, xnew):
    print(f'max({xo}, {xn}) = {max(xo, xn)}')



# print(xodd1)
# print(xodd2)
