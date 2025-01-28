project(
    'game_engine',
    'cpp',
    version : '0.0.1',
    default_options : [
#        'buildtype=debugoptimized',
        'buildtype=release',
        'warning_level=3',
        'cpp_std=c++20',
#        'b_sanitize=address'
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
#    'main.cpp',
    'main_single.cpp'
#    'core/memory.cpp',
]

#Intrinics 1: 13866 ns (1.38674 ns/op)
#Intrinics 2: 13837 ns (1.38384 ns/op)
#Intrinics 3: 10559 ns (1.05601 ns/op)
#Intrinics 4: 13856 ns (1.38574 ns/op)
#Intrinics 5: 10880 ns (1.08811 ns/op)
#Intrinics 6: 10891 ns (1.08921 ns/op)
#    Regular: 9197 ns (0.919792 ns/op)
#Intrinics 1: 13886 ns (1.38874 ns/op)
#Intrinics 2: 13826 ns (1.38274 ns/op)
#Intrinics 3: 10550 ns (1.05511 ns/op)
#Intrinics 4: 10490 ns (1.0491 ns/op)
#Intrinics 5: 10830 ns (1.08311 ns/op)
#Intrinics 6: 10871 ns (1.08721 ns/op)
#    Regular: 16481 ns (1.64826 ns/op)

#thread_dep = dependency('OpenMP')

#Intrinics 1: 55925 ns (0.559256 ns/op)
#Intrinics 2: 55485 ns (0.554856 ns/op)
#Intrinics 3: 49584 ns (0.495845 ns/op)
#Intrinics 4: 48281 ns (0.482815 ns/op)
#Intrinics 5: 49283 ns (0.492835 ns/op)
#Intrinics 6: 48751 ns (0.487515 ns/op)
#    Regular: 82576 ns (0.825768 ns/op)


game_engine = executable(
    'game_engine',
    sources: sources,
    cpp_args: gcc_opts,
    include_directories: includes,
#    dependencies: thread_dep
)

#test('test', game_engine)