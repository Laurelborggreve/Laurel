# -*- coding: utf-8 -*-
"""
DCC_sim.py

Purpose:
    Purpose = Simulate data from a multivariate GARCH DCC model

    Model = The Dynamic Conditional Correlation (DCC) model is given by:
    H(t) = D(t) * R(t) * D(t)
    D(t)^2 = diag{w(i)} + diag{k(i)} * (r(t-1)*r(t-1)') * (diag{labda(i)}*D(t-1)^2)     OR    D(t) = diag{sqrt(h(it))}
    Q(t) = (1 - alpha - beta) * S + (alpha * eps(t-1) * eps(t-1)') + (beta * Q(t-1))
    R(t) = diag{Q(t)}^-0.5 * Q(t) * diag{Q(t)}^-0.5


Version:
    0

Date:
    28-04-2018

@author: Laurel Borggreve
Msc Econometrics and Operations Research
Thesis: Active vs Passive Investing
"""
###########################################################
### Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as pltasdasdasd

###########################################################
### DCC(dOmega, dAlpha, dBeta, iN, mS):

def DCC(dOmega, dAlpha, dBeta, iN, mS):
    """
    Purpose:
      Simulate

    Inputs:
      ...

    Return value:
      mH      conditional covariance matrix
      mR      conditional correlation matrix of standardized disturbances eps(t)
      mS      unconditional covariance of epsilons

    """
     mH = np.zeros(shape=(iN,iN))
     mQ = np.zeros(shape=(iN,iN))
     mR = np.zeros(shape=(iN,iN))


    # Estimation
    for t in range(iN):
        mQ[t] = (1 - dAlpha - dBeta) * mS + (dAlpha * (vEps[t-1] * np.transpose(vEps[t-1]))) + (dBeta * mQ[t-1])
        mR[t] = np.diagonal((1/mQ[t]) * mQ[t] * np.diagonal(1/mQ[t]))
        mH[t] = mD[t] * mR[t] * mD[t]

    return (mQ, mR, mH)


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

def mS(iN, vS2):
    # Create mD
    mD = np.zeros((iN,iN), int)

    for t in range(iN):
        mD[t] = np.fill_diagonal(mD, vS2)

    return (mD)

    # Create vEps
    vEps = vS2 * vY

    # Create mS
    mS = np.zeros(shape=(iN,iN))
    mS = np.cov(vEps, (np.transpose(vEps)))
    return mS

###########################################################
### Output(vY, vS2)
 #def Output(vY, mH):
    # Plot data
  #  fig = plt.figure()

   # plt.subplot(2,2,1)
   # plt.plot(vY, label="returns")
   # plt.plot(2*np.sqrt(vS2), label="2sd")
    # x = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000] # I am trying to maken the plots larger and change the x-as values, but it doesn't work
    # plt.xlabel('time')
    # plt.ylabel('returns')
   # plt.legend()
   # plt.title('Plot data')

    # plt.subplot(2,2,3)
    # plt.plot(np.sqrt(vSig[1:iN]))
    # x = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    # plt.ylabel('y-axis')
    # plt.subplots_adjust(left=0.8, bottom=0.8, right=0.9, top=0.9, wspace=0.8, hspace=0.9)

  #  plt.show()


###########################################################
### main
def main():
    # Magic numbers
    iN = 10000
    dOmega = 0.1
    dAlpha = 0.05
    dBeta = 0.94

    # GenrGARCH
    (vY, vS2)= GenrGARCH(dOmega, dAlpha, dBeta, iN)

    # compute mS
    mS = mS(iN, vS2, vY)

    # Estimation - DCC
    (mQ, mR, mH)= dCC(dOmega, dAlpha, dBeta, iN, mS)

    # Output
    Output(mQ, mR, mH)

###########################################################
### start main
if __name__ == "__main__":
    main()
