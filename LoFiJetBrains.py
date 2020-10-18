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
import numpy as np


class Model:
    """This Class Contains all Methods Relevant to the Multi Fidelity Analysis of the Instantiated Model"""

    # Constructor method
    # ***************************************************************************************************
    def __init__(self, length: float, base: float, height: float, modulus: float):

        self.length = length
        self.base = base
        self.height = height
        self.modulus = modulus

    # ***************************************************************************************************

    # Class-defined methods
    # ***************************************************************************************************
    # average of different experimental data sets under the same load
    def dataAverage(self, data1, data2, data3, data4, data5):
        averageData = np.zeros(len(data1))
        for x in range(len(data1)):
            averageData[x] = (data1[x] + data2[x] + data3[x] + data4[x] + data5[x])/5
        return averageData

    # Calcualtes moment of inertia of Cantilever Beam
    # Enter values in meters
    def inertiaCalc(self, base, height):
        inertia = (1 / 12) * base * (height ** 3)
        return inertia

    # calculates theoretical deflection vs. load
    def theoDataCalc(self, length, inertia, modulus, load):
        count = len(load)
        theoData = []
        for x in range(count):
            entry = ((length ** 4) / (8 * modulus * inertia) * load[x])
            theoData.append(entry)

        return theoData

    # quantifies uncertainty in simulated data with respect to high fidelity experimental data
    def uncertaintyCal(self, avgExpData, simData):
        count = len(simData)
        error = []
        exp = np.array(avgExpData)
        sim = np.array(simData)

        for x in range(count):
            val = exp[x] - sim[x]
            # Magnitude of error
            error.append(val)

        return error

    # approximate error as a linear function varying with load
    # returns 2x1 matrix of coefficients, the first being constant and the second the load coefficient
    def leastSquaresCalc(self, error, load):
        e = np.array(error)
        w = np.array(load)
        count = len(e)
        k = (count, 2)
        A = np.zeros(k)

        for x in range(count):
            A[x, 0] = 1
            A[x, 1] = w[x]

        AT = np.transpose(A)

        def multiply(matrix1, matrix2):
            s = (len(matrix1), len(matrix2[0]))
            res = np.zeros(s)
            for i in range(len(matrix1)):
                for j in range(len(matrix2[0])):
                    res[i, j] = np.dot(matrix1[i, :], matrix2[:, j])
            return res

        ATA = multiply(AT, A)
        ATAinv = np.linalg.inv(ATA)
        ATAinvAT = multiply(ATAinv, AT)
        C1 = np.dot(ATAinvAT[0, :], e)
        C2 = np.dot(ATAinvAT[1, :], e)
        C = [C1, C2]
        return C
        # ***************************************************************************************************
