from scripts.pyforce import *

smith_data = '/home/cepheid/TriForce/game_engine/tests/smith_data'

data_paths = [
    '/home/cepheid/TriForce/game_engine/data/lsi_og',
    '/home/cepheid/TriForce/game_engine/data/rlsi_og',
    '/home/cepheid/TriForce/game_engine/data/rlsi_og_10eV',
    '/home/cepheid/TriForce/game_engine/data/rlsi_og_rad',
    '/home/cepheid/TriForce/game_engine/data/rlsi_og_rad_2ndorder',
    '/home/cepheid/TriForce/game_engine/data/rlsi_og_rc',
    '/home/cepheid/TriForce/game_engine/data/rlsi_og_rc_2ndorder',
    # '/home/cepheid/TriForce/game_engine/data/rlsi_dt',
    # '/home/cepheid/TriForce/game_engine/data/rlsi_dt_rad',
    # '/home/cepheid/TriForce/game_engine/data/rlsi_og_rad_2ndorder',
    # '/home/cepheid/TriForce/game_engine/data/rlsi_dt_rc',
    # '/home/cepheid/TriForce/game_engine/data/rlsi_og_rc_2ndorder',
    # '/home/cepheid/TriForce/game_engine/data/rlsi_dt_rcf',
]

plot_energies(data_paths, smith_data)
plot_spectra(data_paths, smith_data)

# groups = ['deuterium', 'electrons', 'photons', 'neutrons', 'protons', 'helium3', 'tritium']
# for p in data_paths:
#     plot_density(groups, 15000, p, block=True)















