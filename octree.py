#!/usr/bin/env python3

import numpy as np
from functools import cmp_to_key

nx = 7
ny = 1
nz = 6

def mortonEncode(x, y, z):
    answer = 0
    for ii in range(21):
        answer |= ((x & (1 << ii)) << 2 * ii) | ((y & (1 << ii)) << (2 * ii + 1)) | ((z & (1 << ii)) << (2 * ii + 2))
    return answer

codes = []

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            codes.append(mortonEncode(k, j, i))

codes = np.asarray(codes)

print(codes.reshape((nx, ny, nz)))


