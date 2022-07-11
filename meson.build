project('LATfield2','cpp','c',
        default_options : ['cpp_std=c++17',
        'warning_level=3'],
        version: '1.1')

cpp=meson.get_compiler('cpp')
mpi=dependency('mpi',language: 'cpp')
boost = dependency('boost', modules: ['mpi','serialization','unit_test_framework'],required: false)
mpirun = find_program('mpirun')
pkg=import('pkgconfig')



conf_data = configuration_data({ 
    'FFT3D' :get_option('FFT3D'),
    'HDF5' :get_option('HDF5'),
    'SINGLE' :get_option('SINGLE'),
    'EXTERNAL_IO' :get_option('EXTERNAL_IO')
    })

deps = [mpi,boost]
fftw3 = dependency('fftw3', required: get_option('FFT3D'))
if (fftw3.found() and get_option('FFT3D'))
    deps += fftw3
endif
hdf5 = dependency('hdf5', required: get_option('HDF5'))
if (hdf5.found() and get_option('HDF5'))
    deps += hdf5
endif

subdir('include')
subdir('src')

liblatfield = library('latfield',
    latfield_headers,
    part_headers,
    latfield_hdr,
    sources,
    config,
    include_directories: include,
    install: true,
    dependencies: deps)

install_headers(latfield_hdr)
install_headers(latfield_headers, subdir: 'LATfield2')
install_headers(config, subdir: 'LATfield2')
install_headers(part_headers,subdir: 'LATfield2/particles')

pkg.generate(liblatfield,
    description: 'LATfield2')

subdir('benchmarks_tests')
subdir('examples')

subdir('tests')
