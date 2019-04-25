"""
Program: Dynamical Systems Template
File: dynamicalSystemsTemplate.pylab
Author: C.D. Wentworth
Version: 2-10-2018.1
Summary:
    Basic script for solving a system of first order differential equations
    that describes a dynamical system.
Usage: python dynamicalSystemsTemplate.py
Version History:
    1-15-2018.1: base
    12-11-2018.1: added some documentation
    
"""

import scipy.integrate as si
import numpy as np
import matplotlib.pylab as plt

# Create a function that defines the rhs of the differential equation system
def calc_RHS(x,t):
   #y = a list that contains the system state
   #t = the time for which the right-hand-side of the system equations
    #   is to be calculated.
   #p = a tuple that contains any parameters needed for the model
    
    a = 58.0
    b = 7.9e-5
    c = 125.0
    du = 0.04
    di = 2.0
    dv = 3.1
    
   #Unpack the state of the system
    Tu = x[0]
    Ti = x[1]
    V = x[2]

#   Calculate the rates of change (the derivatives)
    dTudt = a - b*Tu*V - du*Tu
    dTidt = b*Tu*V - di*Ti
    dVdt = c*Ti - dv*V

    return [dTudt,dTidt,dVdt]

HIVData = np.loadtxt('HIVData.txt',skiprows = 2)
Week = HIVData[:,0]
TuData = HIVData[:,1]
VData = HIVData[:,2]

# Define the initial conditions
x_0 = (1100,0,10)
# Define the time grid
t = np.linspace(0,50,200)

# Solve the DE
sol = si.odeint(calc_RHS,x_0,t)
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
