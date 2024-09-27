project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
#        'buildtype=debug',
        'warning_level=3',
        'cpp_std=c++20'
    ]
)

gcc_opts = [
# '-fconcepts-diagnostics-depth=200',
 '-march=native',
 '-Wconversion',
# '-save-temps'
 ]

includes = include_directories('core', 'electromagnetics', 'aydenstuff')

sources = [
    'main.cpp',
#    'core/memory.cpp',
]

#thread_dep = dependency('OpenMP')

game_engine = executable(
    'game_engine',
    sources: sources,
    cpp_args: gcc_opts,
    include_directories: includes,
#    dependencies: thread_dep
)

#test('test', game_engine)