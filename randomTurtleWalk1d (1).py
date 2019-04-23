"""

Title: Random Turtle Walk in 1D
Date: 12/30/2018
Author: C.D. Wentworth
version: 12.30.2018.1

"""
import numpy as np
import random as rn
import turtle as trt

def step(xi,lx,ux):
    import random as rn
    r = rn.random()
    if r < 0.5:
        if xi > lx:
            xi = xi -1
    else:
        if xi < ux:
            xi = xi+1
    return xi

#--Main Program
# set up the turtle window
wn = trt.Screen()
lx = -10
ux = 10
ly = -5
uy = 5
wn.setworldcoordinates(lx,ly,ux,uy) 

# create a turtle
t = trt.Turtle()
t.speed(3)
xi = 0
yi = 0

# set up random generator
rn.seed(413)
N = 200
# execute random walk
xlist = [xi]
for i in range(N):
    xi=step(xi,lx,ux)
    xlist.append(xi)
    t.goto(xi,yi)
xnp = np.array(xlist)
xmean = xnp.mean()
xsd = np.sqrt(xnp.var()/float(N-1))

print('mean x = ',xmean,' sd = ',xsd)

wn.exitonclick() 

