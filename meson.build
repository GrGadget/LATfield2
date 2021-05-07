project('LATfield2','cpp','c',
        default_options : ['cpp_std=c++17',
        'warning_level=3'],
        version: '1.1')
add_global_arguments(['-Ofast','-march=native' ],language : 'cpp')

cpp=meson.get_compiler('cpp')
mpi=dependency('mpi',language: 'cpp')
boost = dependency('boost', modules: ['mpi','serialization','unit_test_framework'],required: false)
mpirun = find_program('mpirun')
pkg=import('pkgconfig')

deps = [mpi]
fftw3 = dependency('fftw3', required: get_option('FFT3D'))
if (fftw3.found() and get_option('FFT3D'))
    deps += fftw3
endif
hdf5 = dependency('hdf5', required: get_option('HDF5'))
if (hdf5.found() and get_option('HDF5'))
    deps += hdf5
endif

include=include_directories('.')

latfield_headers = files(
    'Imag.hpp',
    'int2string.hpp',
    'LATfield2_Catalyst.hpp',
    'LATfield2_macros.hpp',
    'LATfield2_Field.hpp',
    'LATfield2.hpp',
    'LATfield2_IO_server.hpp',
    'LATfield2_Lattice.hpp',
    'LATfield2_parallel2d.hpp',
    'LATfield2_PlanFFT.hpp',
    'LATfield2_save_hdf5.hpp',
    'LATfield2_save_hdf5_pixie.h',
    'LATfield2_SettingsFile.hpp',
    'LATfield2_Site.hpp',
    'timer.hpp')

subdir('particles')
subdir('src')

conf_data = configuration_data({ 
    'FFT3D' :get_option('FFT3D'),
    'HDF5' :get_option('HDF5'),
    'H5_HAVE_PARALLEL' :get_option('H5_HAVE_PARALLEL'),
    'SINGLE' :get_option('SINGLE'),
    'EXTERNAL_IO' :get_option('EXTERNAL_IO')
    })
config=configure_file(
    input: 'config.h.in', 
    output: 'config.h',
    configuration: conf_data)

liblatfield = library('latfield',latfield_headers,part_headers,sources,config,
    include_directories: include,
    install: true,
    dependencies: deps)

install_headers(latfield_headers)
install_headers(config)

pkg.generate(liblatfield,
    description: 'LATfield2')

subdir('benchmarks_tests')
subdir('examples')

subdir('tests')
