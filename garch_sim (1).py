#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
garch_sim.py

Purpose:
    Purpose = Simulate data from a univariate GARCH(1,1)

    Model = Generalized Autoregressive Conditional Heteroskedasticity (GARCH(1,1)) model given by:
    x(t) = sqrt(sig(t)) * eps(t)
    sig(t+1) = omega + alpha * x(t)^2 + beta * sig(t)
    where eps is iid distributed > eps ~ NID(0,1)    ...

Version:
    1       with edits from CS Bos

Date:
    21-04-2018

@author: Laurel Borggreve
Msc Econometrics and Operations Research
Thesis: Active vs Passive Investing
"""
###########################################################
### Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

###########################################################
### GenrGARCH(dOmega, dAlpha, dBeta, iN):
def GenrGARCH(dOmega, dAlpha, dBeta, iN):
    """
    Purpose:
      Simulate

    Inputs:
      ...

    Return value:
      vY      iN vector with observations
      vS2     iN vector with variances
    """
    dS2 = dOmega/(1 - dAlpha - dBeta)

    vEps = np.random.normal(size=iN);

    # Define Time Series Vector
    vS2 = np.zeros_like(vEps);
    vY = np.zeros_like(vEps);

    # Estimation
    for t in range(iN):
        vY[t] = np.sqrt(dS2) * vEps[t]
        vS2[t]= dS2
        # Update
        dS2 = dOmega + dAlpha * vY[t]**2 + dBeta * dS2

    return (vY, vS2)


###########################################################
### Output(vY, vS2)
def Output(vY, vS2):
    # Plot data
    fig = plt.figure()

    plt.subplot(2,2,1)
    plt.plot(vY, label="returns")
    plt.plot(2*np.sqrt(vS2), label="2sd")
    # x = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000] # I am trying to maken the plots larger and change the x-as values, but it doesn't work
    # plt.xlabel('time')
    # plt.ylabel('returns')
    plt.legend()
    plt.title('Plot data')

    # plt.subplot(2,2,3)
    # plt.plot(np.sqrt(vSig[1:iN]))
    # x = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    # plt.ylabel('y-axis')
    # plt.subplots_adjust(left=0.8, bottom=0.8, right=0.9, top=0.9, wspace=0.8, hspace=0.9)

    plt.show()


###########################################################
### main
def main():
    # Magic numbers
    iN = 10000
    dOmega = 0.1
    dAlpha = 0.05
    dBeta = 0.94

    # Estimation
    (vY, vS2)= GenrGARCH(dOmega, dAlpha, dBeta, iN)

    # Output
    Output(vY, vS2)

###########################################################
### start main
if __name__ == "__main__":
    main()


###########################################################
### create mD
    







