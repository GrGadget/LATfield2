# mass = executable('mass',['mass.cpp'],
#     include_directories: include,
#     dependencies: [boost,boost_ut])
# test('particle mass with SFINAE',mass)
#     
# if (boost.found())
#     fft = executable('fft',['fft.cpp'],
#         include_directories: include,
#         dependencies: [boost,mpi,hdf5,fftw3],
#         link_with: liblatfield)
#     test('FFT correctedness',mpirun,
#        args: ['-np','4',fft.full_path()],
#        is_parallel: false)
# endif
