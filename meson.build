project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
#        'buildtype=debug',
        'buildtype=debugoptimized',
#        'buildtype=release',
        'warning_level=3',
        'cpp_std=c++23',
#        'b_sanitize=address',
        'b_lto=true',
#        'b_pgo=generate',
#        'b_pgo=use'
    ]
)

debug_opts = [
#    '-O0', '-g'
#    '-fopt-info-vec-missed=missed.txt',
#    '-DNDEBUG',
#    '-fno-omit-frame-pointer',
#    '-fno-inline-functions',
#    '-fno-inline-functions-called-once',
#    '-fno-optimize-sibling-calls'
]

gcc_opts = [
    '-Wconversion',
    '-fuse-ld=mold',
    '-march=native',
    '-ftree-vectorize',
#    '-ffast-math',
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
    'electromagnetics/src/em_solver.cpp',
    'electromagnetics/src/sources.cpp',
    'metrics/metrics.cpp',
    'particles/src/particles.cpp',
    'particles/src/pusher.cpp'
]


deps = [
    dependency('adios2', required: true),
    dependency('OpenMP')
]

game_engine = executable(
    'game_engine',
    sources: sources,
    cpp_args: gcc_opts + debug_opts,
    include_directories: includes,
    dependencies: deps
)
#test('test', game_engine)