import numpy as np

data_path = '/home/cepheid/TriForce/game_engine/tests/mikes_files/DD_pT_BH_keV_mbarn_231009.dat'
output_path = '/tests/collision_tests/cross_section_data/DD_pT_BH_eV_m2.txt'

data = np.genfromtxt(data_path)

data[:, 0] *= 1000.0 # keV to eV
data[:, 1] *= 1.0e-31 # millibarn to m^2

np.savetxt(output_path, data)

