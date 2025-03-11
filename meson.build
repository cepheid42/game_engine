project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
        'buildtype=debugoptimized',
#        'buildtype=release',
        'warning_level=3',
        'cpp_std=c++23',
#        'b_sanitize=address'
    ]
)

gcc_opts = [
# '-fconcepts-diagnostics-depth=200',
 '-march=native',
 '-Wconversion',
# '-save-temps'
 ]

includes = include_directories('core', 'electromagnetics', 'metrics')

sources = [
    'main.cpp',
    'electromagnetics/em_solver.cpp',
    'electromagnetics/em_sources.cpp',
    'metrics/metrics.cpp'
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