# UndariaGEN - A spatially-explicit agent-based model of macroalgae in the marine environment
Copyright (C) 2005-2015  James T. Murphy, Ray Walshe

This is a individual-based model of growth and dynamics of the invasive macroalgal species Undaria pinnatifida. The model includes explicit representations of the various life stages of U. pinnatifida (microscopic gametophytes and macroscopic sporophytes) and their responses to environmental parameters such as light and temperature using empirical data from the literature.  

The model framework can be used to explicitly represent complex spatial and temporal patterns of invasion in order to be able to make quantitative predictions about the impact of these factors on invasion dynamics of U. pinnatifida. This would be a useful tool for making risk assessments of invasion potential under different environmental conditions and for choosing optimal control strategies for cost-effective management.


#1. Compilation

Requires: 
   - MPICH (Implementation of the MPI Message Passing Interface):
      - http://www.mpich.org/
      - Ubuntu: install mpich and libmpich-dev packages
   - GNU Scientific Library (GSL):
      - http://www.gnu.org/software/gsl/
      - Ubuntu: install gsl-bin and libgsl0-dev packages

To compile:
   - make
   

#2. Execution

To run:
   - mpiexec -n 4 ./bin/undariagen [./input/inFile1.in] [./input/inFile2.in]

Example run:
   - mpiexec -n 4 ./bin/undariagen ./input/test_brest.in ./input/port_edit514_482.in

Input files:
   - inFile1.in = input parameters (see Sec. 3) [test_brest.in]
   - inFile2.in = map of substrate types (imageJ xy coords) [port_edit514_482.in]
