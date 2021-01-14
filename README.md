# Scalar-Field-Monte-Carlo

## Introduction
This code is designed to run simple Monte Carlo simulations for Scalar Field Theory simulations.
It numerically computes observables through the Path Integral and outputs the observables to a text file.
Tools for analyzing the data are included in a Mathematica Notebook.  If you choose not to use these,
know that the code produces data for every field configuration after the first 1000 and that all these numbers
will be highly correlated.  Determining the autocorrelation and then either spacing out samples or using binning 
is the wisest decision to deal with these issues.


## Purpose
This code was written as an introductory research project and as a means to study methods of symmetry breaking in
Field theories.  It is currently written to output the average field configuration of a complex scalar field, along
with the action at each field point.  The ultimate goal is to study ways to overcome the sign problem for finite density
simulations.  But for now, it serves as my way of learning more about symmetries and coding.


## Acknowledgements
This code uses the Mersenne Twister random number generator.  More information can be found in the mt64 source file.  
Additionally, The data analysis routines were given to me by Dr. Andrei Alexandru of The George Washington University.

