latfield_headers = files(
    'Imag.hpp',
    'int2string.hpp',
    'Catalyst.hpp',
    'macros.hpp',
    'Field.hpp',
    'IO_server.hpp',
    'Lattice.hpp',
    'parallel2d.hpp',
    'PlanFFT.hpp',
    'save_hdf5.hpp',
    'save_hdf5_pixie.hpp',
    'SettingsFile.hpp',
    'Site.hpp',
    'timer.hpp')
subdir('particles')
config=configure_file(
    input: 'config.h.in', 
    output: 'config.h',
    configuration: conf_data)
