# mass = executable('mass',['mass.cpp'],
#     include_directories: include,
#     dependencies: [boost,boost_ut])
# test('particle mass with SFINAE',mass)
#     
fft = executable('fft',['fft.cpp'],
    include_directories: include,
    dependencies: [boost,boost_ut,mpi,hdf5,fftw3,boost_mpi],
    link_with: liblatfield)
test('FFT correctedness',mpirun,
   args: ['-np','4',fft.full_path()],
   is_parallel: false)
