gpusvcalibration
================

R package to support fast calibration of stochastic volatility models for option pricing using GPUs

Installation instructions:

Download the linux x64 binary and install it with the command:

install.packages('gpusvcalibration_0.0-1.tar.gz', repo=NULL)

or clone the repository and then

Modify the R_EXE and R_INCLUDE variables in the Makefile

R CMD build gpusvcalibration

R CMD check gpusvcalibration

R CMD install gpusvcalibration_0.0-1.tar.gz <target-filepath>

install.packages('<target-filepath>gpusvcalibration_0.0-1.tar.gz', repo=NULL)
