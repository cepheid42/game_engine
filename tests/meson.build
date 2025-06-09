
# Single particle moving in a constant +-1 V/m E-field
sp_accel_test = executable(
    'sp_accel_test',
    sources: ['/home/cepheid/TriForce/game_engine/tests/sp_accel_test.cpp'],
    cpp_args : gcc_opts + debug_opts,
    include_directories : includes,
    dependencies : deps
)

# Single electron moving in a constant 4.2e-4 T B-field
sp_gyro_test = executable(
    'sp_gyro_test',
    sources: ['/home/cepheid/TriForce/game_engine/tests/sp_gyro_test.cpp'],
    cpp_args : gcc_opts + debug_opts,
    include_directories : includes,
    dependencies : deps
)

# Tests single particle deposition
sp_gyro_test = executable(
    'sp_jdep_test',
    sources: ['/home/cepheid/TriForce/game_engine/tests/sp_jdep_test.cpp'],
    cpp_args : gcc_opts + debug_opts,
    include_directories : includes,
    dependencies : deps
)