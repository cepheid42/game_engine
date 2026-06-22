import numpy as np

data_path = '/home/cepheid/TriForce/game_engine/tests/mikes_files/SB_G4_Z1_kdsdk_MeV_barns.csv'
output_path = '/tests/collision_tests/cross_section_data/SB_G4_Z1_kdsdk_eV_m2.txt'

data = np.genfromtxt(data_path)

data[:, 0] *= 1000.0 # keV to eV
data[:, 1] *= 1.0e-31 # millibarn to m^2

np.savetxt(output_path, data)

