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

# dataCorrection(error: list, controlVar: list, oldData: list, order:int) = errorModel, newData, covariance

# *******************************************************************************************************
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math

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

    # Calculates moment of inertia of Cantilever Beam
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
    def uncertaintyCal(self, LoFiData, HiFiData):
        count = len(HiFiData)
        error = np.zeros(count)
        LoFi = np.array(LoFiData)
        HiFi = np.array(HiFiData)

        for i in range(count):
            val = np.abs(LoFi[i] - HiFi[i]) /HiFi[i]
            # Magnitude of error
            error[i] = val
        return error

    def dataCorrection(self, error, controlVar, oldData, order):
        count = len(error)

        if order == 1:

            def linear(x, a, b):
                return(a*x + b)

            poptLinear, pcovLinear = curve_fit(linear, controlVar, error)
            a = poptLinear[0]
            b = poptLinear[1]

            newData = np.zeros(count)
            calcError = np.zeros(count)
            for i in range(count):
                calcError[i] = linear(controlVar[i], a, b)
                newData[i] = oldData[i] - calcError[i]* oldData[i]
            return calcError, newData, pcovLinear

        elif order == 2:
            def quadratic(x, a, b, c):
                return (a * (x ** 2)) + (b * x) + c

            poptQuadratic, pcovQuadratic = curve_fit(quadratic, controlVar, error)
            a = poptQuadratic[0]
            b = poptQuadratic[1]
            c = poptQuadratic[2]

            newData = np.zeros(count)
            calcError = np.zeros(count)
            for i in range(count):
                calcError[i] = quadratic(controlVar[i], a, b, c)
                newData[i] = oldData[i] - calcError[i] * oldData[i]
            return calcError, newData, pcovQuadratic

        elif order == 3:
            def cubic(x, a, b, c, d):
                return (a * (x ** 3)) + (b * (x ** 2)) + (c * x) + d

            poptCubic, pcovCubic = curve_fit(cubic, controlVar, error)
            a = poptCubic[0]
            b = poptCubic[1]
            c = poptCubic[2]
            d = poptCubic[3]
            print(a,b,c,d)

            newData = np.zeros(count)
            calcError = np.zeros(count)
            for i in range(count):
                calcError[i] = cubic(controlVar[i], a, b, c, d)
                newData[i] = oldData[i] - calcError[i]* oldData[i]
            return calcError, newData, pcovCubic

        elif order == 4:
            def quartic(x, a, b, c, d, e):
                return (a * (x ** 4)) + (b * (x ** 3)) + (c * (x ** 2)) + (d*x) + e

            poptQuartic, pcovQuartic = curve_fit(quartic, controlVar, error)
            a = poptQuartic[0]
            b = poptQuartic[1]
            c = poptQuartic[2]
            d = poptQuartic[3]
            e = poptQuartic[4]

            newData = np.zeros(count)
            calcError = np.zeros(count)
            for i in range(count):
                calcError[i] = quartic(controlVar[i], a, b,c ,d ,e)
                newData[i] = oldData[i] - calcError[i]* oldData[i]
            return calcError, newData, pcovQuartic

        elif order == 5:
            def log(x, a, b):
                return a*np.log(x) + b

            poptLogrithmic, pcovLogrithmic = curve_fit(log, controlVar, error)
            a = poptLogrithmic[0]
            b = poptLogrithmic[1]
            print(a,b)

            newData = np.zeros(count)
            calcError = np.zeros(count)
            for i in range(count):
                calcError[i] = log(controlVar[i], a, b)
                newData[i] = oldData[i] - calcError[i]* oldData[i]
            return calcError, newData, pcovLogrithmic

        elif order == 6:
            def exp(x, a, b, c):
                return a*(np.e ** (c*x)) + b

            poptExp, pcovExp = curve_fit(exp, controlVar, error)
            a = poptExp[0]
            b = poptExp[1]
            c = poptExp[2]

            newData = np.zeros(count)
            calcError = np.zeros(count)
            for i in range(count):
                calcError[i] = exp(controlVar[i], a, b, c)
                newData[i] = oldData[i] - calcError[i]* oldData[i]
            return calcError, newData, pcovExp

        else:
            print("invalid function")

# ***************************************************************************************************
