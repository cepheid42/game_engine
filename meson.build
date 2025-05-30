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
#        '-O0', '-g',
#        '-fopt-info-note-vec-missed=missed.txt',
#        '-DNDEBUG',
#        '-fno-omit-frame-pointer',
#        '-fno-inline-functions',
#        '-fno-inline-functions-called-once',
#        '-fno-optimize-sibling-calls',
#        '-fno-default-inline',
#        '-fno-inline'
]

gcc_opts = [
    '-Wconversion',
    '-fuse-ld=mold',
    '-march=native',
    '-ftree-vectorize',
]

includes = include_directories(
    'electromagnetics',
    'metrics',
    'params',
    'particles',
    'utilities'
)

sources = [
    'main.cpp'
]

deps = [
    dependency('adios2', required : true),
    dependency('OpenMP'),
    #    dependency('cuda')
]

game_engine = executable(
    'game_engine',
    sources : sources,
    cpp_args : gcc_opts + debug_opts,
    include_directories : includes,
    dependencies : deps
)

single_test = executable(
    'single_particle',
    sources: ['./tests/single_particle.cpp'],
    cpp_args : gcc_opts + debug_opts,
    include_directories : includes,
    dependencies : deps
)
