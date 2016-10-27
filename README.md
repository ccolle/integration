
WARNING, FIRST VERSION, so do not assume main interface of integrator will not change in future versions...

A start of C++ collection of integration libraries. All wrapped to same common interface
declared in inc/integrator.hpp

Hopefully this makes your Code much more flexible towards changing integration methods and
settings (as precision, max number of evals, etc...)

If you wish to add your own integration routines, see "*_wrapper.hpp" in inc/

[IMPORTANT]
    The Cuba routines are intended for high dimensional integrands!
    For dimensions 1,2,(3?) you are probably better of using the other
    algorithms... (Vegas might work though)

[DEBUGGING]
-Cuba v4.2:
    It seems that valgrind finds errors for the supplied demo-c.c file
    for Suave and Divonne. Vegas and Cuhre seem to work fine though.
    The valgrind warning '72,704 bytes in 1 blocks are still reachable in loss record'
    is apparently a well known bug (or feature according to the devs...)
    
    Cuhre seems to fail miserably for 1D integrands. It does not evaluate
    the integrand at hand once!
    
    


You can always ask Camille for some help...


To build the code you need CMake.

  1-  $> cd build;

  2-  $> cmake -DCMAKE_BUILD_TYPE=RELEASE ..

  3-  $> make;

  4-  $> make install;


Now you are all set. To use this in your code, include the headers in inc/
and link to the libraries in lib/

Your compilation command might look something like (untested...)

g++ yourfile.cpp -I/path/to/integration/inc -L/path/to/integration/lib -Wl,-rpath=/path/to/integration/lib -lcuba -lcubature $(gsl-config --libs) -o yourexec


