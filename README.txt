Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
See the LICENSE.txt file at the top-level of this distribution.

THIS IS BETA SOFTWARE. FEATURES and FUNCTIONALITY ARE SUBJECT TO CHANGE WITHOUT
PRIOR NOTICE.

Description
===========
Library implementation of a probability flux-based method for extracting
momentum distributions from traveling wavepackets in time-dependent quantum
mechanical calculations. Both a semi-classical method and a method
incorporating first-order quantum corrections are implemented.

Required dependencies
=====================
Build requirements:

* Fortran compiler (most recently tested on GNU Fortran 6.2)
* SCons (most recently tested on SCons 2.5)

External libraries to build from source (build scripts included):

* Fortran generics library, available at:
  https://github.com/kramer314/fortran-lib
  Note that this library is still in beta status, so the API is technically
  still in flux. Always use the most recent version.

Recommended dependencies
========================
* OpenMP implementation (most recently tested with GCC OpemMP 6.2)


