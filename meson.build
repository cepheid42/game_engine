project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
        'cpp_std=c++23',
        'warning_level=3',
#        'optimization=0',
#        'debug=true',
        'b_lto=true',
        'optimization=3',
        'debug=false',
#        'werror=true',
#        'b_sanitize=address',
#        'b_pgo=generate',
#        'b_pgo=use'
    ]
)

debug_opts = [
    '-Wconversion',
#    '-Wnrvo',
#    '-fopt-info-note-vec-missed=missed.txt',
#    '-DNDEBUG',
#    '-fno-omit-frame-pointer',
#    '-fno-inline-functions',
#    '-fno-inline-functions-called-once',
#    '-fno-optimize-sibling-calls',
#    '-fno-default-inline',
#    '-fno-inline'
]

gcc_opts = [
#    '-fuse-ld=mold',
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
    'main.cpp',
]

deps = [
    dependency('adios2', required : true),
    dependency('OpenMP', required : true),
#    dependency('cuda')
]

game_engine = executable(
    'game_engine',
    sources : sources,
    cpp_args : gcc_opts + debug_opts,
    include_directories : includes,
    dependencies : deps
)

#subdir('tests')
