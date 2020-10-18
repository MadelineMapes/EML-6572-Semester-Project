# EML-6572-Semester-Project
# Created by Madeline Mapes at the University of Central Florida
# Revised 18 Oct 2020

# LoFiJetBrains is a toolkit for improving predictions made by finite element analysis software by conducting a
# multi fidelity analysis that combines data from simulated results with existing experimental results. A least
# squares approximation is done to estimate the error in the simulated data as it varies with load.

# The current iteration of this toolkit is applicable only to simple, steel, cantilever beams under a distributed load,
# although future iterations will have more diverse applicability with the results of this iteration serving as a basis
# for future analysis.

# requires numpy

# all units are mks

# *******************************************************************************************************
# Methods, inputs, and return types:

# dataAverage(experimentalData1: list, ..., experimentalData5: list) = averageData: list

# inertiaCalc(length: float, base: float, height: float) = inertia: float

# theoDataCalc(length: float, inertia: float, youngsModulus: float, load: list) = theoData: list

# uncertaintyCalc(averageExperimentalData: list, simulatedData: list) = error: list

# leastSquaresCalc(error: list, load: list) = leastSquaresCoefficients: 2x1 matrix

# *******************************************************************************************************
