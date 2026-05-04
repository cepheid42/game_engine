project(
    'tflink3',
    'cpp',
    version : '0.0.1',
    default_options : [
        'optimization=2',
        'debug=true',
        'warning_level=3',
#        'werror=true',
        'cpp_std=c++23',
#        'b_sanitize=address',
        'b_lto=true',
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
    dependency('OpenMP'),
    dependency('openPMD'),
]

executable(
    'tflink3',
    sources : sources,
    cpp_args : gcc_opts + debug_opts,
    include_directories : includes,
    dependencies : deps
)

#subdir('tests')
