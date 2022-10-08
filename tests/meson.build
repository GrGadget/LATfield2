# mass = executable('mass',['mass.cpp'],
#     include_directories: include,
#     dependencies: [boost,boost_ut])
# test('particle mass with SFINAE',mass)
#     
if (boost.found())
    fft = executable('fft',['fft.cpp'],
        include_directories: include,
        dependencies: [boost,mpi,hdf5,fftw3],
        link_with: liblatfield)
    test('FFT correctedness',mpirun,
       args: ['--oversubscribe','-np','4',fft.full_path()],
       is_parallel: false)
    
    cart_comm = executable('cart_comm',['cart_comm.cpp'],
        include_directories: include,
        dependencies: [boost,mpi,hdf5,fftw3],
        link_with: liblatfield)
    test('Cartesian communicator',mpirun,
       args: ['--oversubscribe','-np','12',cart_comm.full_path()],
       is_parallel: false)
endif
