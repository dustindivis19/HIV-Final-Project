"""
Program: HIV Model 2 - Curve Fit
File: HIVModel2_fit.py
Author: C.D. Wentworth
Version: 4.22.2019.1
Summary:
    This program implements a dynamical systems model of the time development
    of the HIV infection in an individual.  It is a 3 compartment model: 
    uninfected T cells, infected T cells, and virus particles. It attempts
    to fit one of the model paraments to the data using the least squares
    technique.
Usage: python HIVModel2_fit.py
Version History:
    4.22.2019.1: base

    
"""

import scipy.integrate as si
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit

# Create a function that defines the rhs of the differential equation system
def calc_RHS(x,t,p):
   #y = a list that contains the system state
   #t = the time for which the right-hand-side of the system equations
    #   is to be calculated.
   #p = a tuple that contains any parameters needed for the model
    
#    a = 38.11
    a = 38.11
    #b = 8.78e-5
    b = 8.78e-5
    #c = 128.17
    c = 128.17
    #du = 0.08
    du = p[0]
    di = 2.28
    dv = 3.47
    
   #Unpack the state of the system
    Tu = x[0]
    Ti = x[1]
    V = x[2]

#   Calculate the rates of change (the derivatives)
    dTudt = a - b*Tu*V - du*Tu
    dTidt = b*Tu*V - di*Ti
    dVdt = c*Ti - dv*V

    return [dTudt,dTidt,dVdt]

def Tuf(t,du):
    import numpy as np
    import scipy.integrate as si
    Tu0 = x_0[0]
    p = (du,)
    TuL = []
    for tt in t:
        if np.abs(tt)<1.e-5:
            TuLt = Tu0
        else:
            ta = np.linspace(0,tt,10)
            # Solve the DE
            sol = si.odeint(calc_RHS,x_0,ta,args=(p,))
            TuLt = sol[-1,0]
        TuL.append(TuLt)
    return TuL

#--Main Program
    
# Read in data
HIVData = np.loadtxt('HIVData.txt',skiprows = 2)
Week = HIVData[:,0]
TuData = HIVData[:,1]
VData = HIVData[:,2]

# Define the initial conditions
x_0 = [1100,0,10]

# Find optimum model parameter
du = 0.08
p = (du,)
popt,pcov = curve_fit(Tuf,Week,TuData,p)
p_stderr = np.sqrt(np.diag(pcov))

# Calculate theory
# Define the time grid
t = np.linspace(0,50,200)
du = popt[0]
p = (du,)
sol = si.odeint(calc_RHS,x_0,t,args=(p,))
TuTheory = sol[:,0]
TiTheory = sol[:,1]
VTheory = sol[:,2]

# Plot the solution
fig, ax1 = plt.subplots()

ax1.set_xlabel('$t \  [Wks]$')
ax1.set_ylabel('$T-Cell \ Concentration \ [Cells/\mu L]$',color='b')
ax1.plot(t, TuTheory , color='b',label='Tu - Theory')
ax1.plot(t,TiTheory, color='k',label='Ti - Theory')
ax1.tick_params(axis='y',color='b')

ax1.plot(Week,TuData,linestyle='',marker='d',markersize=8.0,color='b',
         label='Tu - Data')

ax2 = ax1.twinx()  # initiate a second axes that shares the same x-axis

ax2.set_ylabel('$HIV \ RNA \ Concentration \ [Copies/\mu L]$',color='r')  # we already handled the x-label with ax1
#ax2.plot(Week, VTheory, color='r')
ax2.tick_params(axis='y',color='r')
ax2.plot(t,VTheory,color='r',label='V - Data')

ax2.plot(Week,VData,linestyle='',marker='^',markersize=8.0,color='r',
         label='V - Theory')
ax1.legend(loc='upper right')
ax2.legend(loc=(.7,0.6))
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.savefig('HIVModel2.png')
plt.show()
print(du)
