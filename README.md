gpusvcalibration
================

R package to support fast calibration of stochastic volatility models for option pricing using GPUs

Installation instructions:

Download the linux x64 binary and install it with the command:

install.packages('gpusvcalibration_0.0-1.tar.gz', repo=NULL)

or clone the repository and then
R CMD build gpusvcalibration
R CMD check gpusvcalibration
R CMD install gpusvcalibration <target-filepath>
