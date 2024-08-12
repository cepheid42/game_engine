project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
        'warning_level=3',
        'cpp_std=c++20'
    ]
)

gcc_opts = [
 '-fconcepts-diagnostics-depth=200',
 '-march=native',
 '-Wconversion'
 ]

includes = include_directories('core', 'electromagnetics')

sources = [
    'main.cpp',
    'core/memory.cpp',
]


game_engine = executable('game_engine', sources: sources, cpp_args: gcc_opts, include_directories: includes)

#test('test', game_engine)