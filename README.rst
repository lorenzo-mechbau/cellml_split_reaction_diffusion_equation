========================================
CellML Split Reaction Diffusion Equation
========================================

This example solves the weak form of the following static advection-diffusion equation, 







Building the example
====================

The fortran version of the example can be configured and built with CMake::

  git clone https://github.com/OpenCMISS-Examples/cellml_split_reaction_diffusion_equation
  mkdir cellml_split_reaction_diffusion_equation-build
  cd cellml_split_reaction_diffusion_equation-build
  cmake -DOpenCMISSLibs_DIR=/path/to/opencmisslib/install ../cellml_split_reaction_diffusion_equation
  make

This will create the example executable "cellml_split_reaction_diffusion_equation" in ./src/fortran/ directory.

Running the example
===================

Fortran version::

  cd ./src/fortran/
  ./cellml_split_reaction_diffusion_equation

Verifying the example
=====================

Results can be visualised by running `visualise.cmgui <./src/fortran/visualise.cmgui>`_ with the `Cmgui visualiser <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.

The following figure shows the solutions and various field variables ........... 




The expected results from this example are available in `expected_results <./src/fortran/expected_results>`_ folder.  

Prerequisites
=============

There are no additional input files required for this example as it is self-contained.

License
=======

License applicable to this example is described in `LICENSE <./LICENSE>`_.


