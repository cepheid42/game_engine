project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
        'buildtype=debugoptimized',
        'warning_level=3',
        'cpp_std=c++20'
    ]
)

includes = include_directories('core', 'electromagnetics')

sources = [
    'main.cpp',
    'core/memory.cpp',
]


game_engine = executable('game_engine', sources: sources, include_directories: includes)

#test('test', game_engine)