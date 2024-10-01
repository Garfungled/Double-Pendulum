import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__))) # to fix relative import errors

import numpy
from pendulum import Pendulum
from datetime import datetime

def updateDelta(deltaF, deltaI, delta, x):
    return x.solution + (deltaI * (delta.solution - x.solution)) / deltaF

def partialLyapunov(deltaF, deltaI):
    return numpy.log(abs(deltaF/deltaI))

def lyapunov(theta1, theta2, epsilon):
    dt = 1/100
    totalTime = 25
    time = numpy.arange(0, totalTime, dt)
    averages = 20
    
    lyapunovExps = []
    x = Pendulum(theta1, theta2)
    delta = Pendulum(theta1 + epsilon, theta2)
    deltaI = numpy.linalg.norm(delta.solution - x.solution)
    
    i = 0
    for t in time:
        i += 1
        x.step(t, dt)
        delta.step(t, dt)
        
        if i % averages == 0:
            deltaF = numpy.linalg.norm(delta.solution - x.solution)
            
            lyapunovExps.append(partialLyapunov(deltaF, deltaI))
            delta.solution = updateDelta(deltaF, deltaI, delta, x)
            deltaI = numpy.linalg.norm(delta.solution - x.solution)
    
    return numpy.sum(lyapunovExps)/totalTime

startTime = datetime.now()
print(lyapunov(1.0, 1.0, 0.001))
print(datetime.now() - startTime)