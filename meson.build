project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
#        'buildtype=debug',
        'buildtype=debugoptimized',
#        'buildtype=release',
        'warning_level=3',
        'werror=true',
        'cpp_std=c++23',
#        'b_sanitize=address'
    ]
)

gcc_opts = [
    '-fuse-ld=mold',
    '-march=native',
    '-Wconversion'
]

includes = include_directories(
    'electromagnetics/include',
    'metrics',
    'params',
    'particles/include',
    'utilities'
)

sources = [
    'main.cpp',
    'electromagnetics/src/em_data.cpp',
    'electromagnetics/src/em_solver.cpp',
    'electromagnetics/src/sources.cpp',
    'metrics/metrics.cpp',
    'particles/src/particles.cpp',
    'particles/src/pusher.cpp'
]

thread_dep = dependency('OpenMP')
adios_dep = dependency('adios2', required: true)

game_engine = executable(
    'game_engine',
    sources: sources,
    cpp_args: gcc_opts,
    include_directories: includes,
    dependencies: [thread_dep, adios_dep]
)

#test('test', game_engine)