
WARNING, FIRST VERSION, so do not assume main interface of integrator will not change in future versions...

A start of C++ collection of integration libraries. All wrapped to same common interface
declared in inc/integrator.hpp

Hopefully this makes your Code much more flexible towards changing integration methods and
settings (as precision, max number of evals, etc...)

If you wish to add your own integration routines, see "*_wrapper.hpp" in inc/

You can always ask camille for some help...


To build the code you need CMake.

$> cd build;
$> make
$> make install


Now you are all set. To use this in your code, include the headers in inc/
and link to the libraries in lib/


