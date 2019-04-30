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
from scipy.optimize import curve_fit

# Create a function that defines the rhs of the differential equation system
def calc_RHS(x,t):
   #x = a list that contains the system state
   #t = the time for which the right-hand-side of the system equation
    
    a = 38.11
    b = 8.78e-5
    c = 128.17
    du = 0.053
    di = 2.29
    dv = 3.47
    
   #Unpack the state of the system
    Tu = x[0]
    Ti = x[1]
    V = x[2]

   #Calculate the rates of change (the derivatives)
    dTudt = a - b*Tu*V - du*Tu
    dTidt = b*Tu*V - di*Ti
    dVdt = c*Ti - dv*V

    return [dTudt,dTidt,dVdt]

def inoc(t,M,NO,r):
    conc = M*N0*np.exp(r*t)/(M + (np.exp(r*t)-1))
    return conc

HIVData = np.loadtxt('HIVData.txt',skiprows = 2)
Week = HIVData[:,0]
TuData = HIVData[:,1]
VData = HIVData[:,2]

HIVTData = np.loadtxt('HIVTreatmentData.txt',skiprows=1)
Days = HIVTData[:,0]
VirusData = HIVTData[:,1]

# Define initial guess for model parameters
M = 1.0
N0 = 10.0
r = 0.020
p = M,N0,r

# Fit data to exponential growth model
popt,pcov = curve_fit(inoc,Days,VirusData,p)
p_stderr = np.sqrt(np.diag(pcov))

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
ax2.plot(t,VTheory,color='r',label='V - Theory')

ax2.plot(Week,VData,linestyle='',marker='^',markersize=8.0,color='r',
         label='V - Data')
ax1.legend(loc='upper right')
ax2.legend(loc=(.7,0.6))
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.savefig('HIVModel2.png')
plt.title('Concentration of HIV and T-Cells Over Time',fontsize=14)
plt.show()
#plot treatment/innoculation model
plt.plot(Days,VirusData,linestyle='',marker='d',color='b',label='Data')
plt.title('Concentration of HIV RNA Post Inoculation',fontsize=14)
plt.xlabel('Days Post Inoculation')
plt.ylabel('HIV RNA Concentration [Copies/mL]',color='b')
M,N0,r = popt
dM,dN0,dr = p_stderr
tTheory = np.linspace(0,45,50)
ConcTheory = inoc(tTheory,M,N0,r)
plt.plot(tTheory,ConcTheory,color='r',label='Model')
plt.legend()
plt.show()
print('M = ',M,'N0 = ',N0,'r = ',r)
