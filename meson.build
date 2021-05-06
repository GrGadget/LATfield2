project('LATfield2','cpp','c',
        default_options : ['cpp_std=c++17',
        'warning_level=3'],
        version: '1.1')
add_global_arguments(['-Ofast','-march=native' ],language : 'cpp')

cpp=meson.get_compiler('cpp')
mpi=dependency('mpi',language: 'cpp')
boost = dependency('boost', modules: ['mpi','serialization','unit_test_framework'],required: false)
mpirun = find_program('mpirun')
fftw3 = dependency('fftw3')
hdf5 = dependency('hdf5')
pkg=import('pkgconfig')

include=include_directories('.')

latfield_headers = files(
    'Imag.hpp',
    'int2string.hpp',
    'LATfield2_Catalyst.hpp',
    'LATfield2_Field.hpp',
    'LATfield2.hpp',
    'LATfield2_IO_server.hpp',
    'LATfield2_Lattice.hpp',
    'LATfield2_parallel2d.hpp',
    'LATfield2_PlanFFT_decl.hpp',
    'LATfield2_PlanFFT.hpp',
    'LATfield2_save_hdf5.h',
    'LATfield2_save_hdf5_pixie.h',
    'LATfield2_SettingsFile.hpp',
    'LATfield2_Site.hpp',
    'looping_macro.hpp',
    'timer.hpp')

subdir('particles')
subdir('src')


liblatfield = library('latfield',latfield_headers,part_headers,sources,
    include_directories: include,
    install: true,
    dependencies: mpi)
install_headers(latfield_headers)
pkg.generate(liblatfield,
    description: 'LATfield2')

subdir('benchmarks_tests')
subdir('examples')

subdir('tests')
