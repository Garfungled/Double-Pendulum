import numpy
from numpy import sin, cos
from scipy.constants import g, pi


# Only doing this for double pendulum, will make n-pendulum later
class Pendulum:
    def __init__(self, theta1, theta2, l1 = 1.0, l2 = 1.0, m1 = 1.0, m2 = 1.0):
        self.m1, self.m2 = m1, m2
        self.l1, self.l2 = l1, l2
        self.solution = numpy.array([theta1, theta2, 0, 0])
        
    def F(self, t, y):
        theta1, theta2, omega1, omega2 = y[0], y[1], y[2], y[3]

        A = numpy.array([[(self.m1 + self.m2)*self.l1, self.m2*self.l2*cos(theta1 - theta2)],[self.l1*cos(theta1 - theta2), self.l2]])
        B = numpy.array([-self.m2*self.l2*omega2*omega2*sin(theta1 - theta2) - (self.m1 + self.m2)*g*sin(theta1), self.l1*omega1*omega1*sin(theta1 - theta2) - g*sin(theta2)])
        alpha = numpy.linalg.inv(A).dot(B)

        return numpy.array([omega1, omega2, alpha[0], alpha[1]])

    def RK4(self, t, y, dt):
        k1 = self.F(t, y)
        k2 = self.F(t + 0.5*dt, y + 0.5*k1*dt)
        k3 = self.F(t + 0.5*dt, y + 0.5*k2*dt)
        k4 = self.F(t + dt, y + dt*k3)
        
        return (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    def position(self):
        theta1, theta2 = self.solution[0], self.solution[1]
        x1, y1 = self.l1*sin(theta1), self.l1*cos(theta1)
        x2, y2 = x1 + self.l2*sin(theta2), y1 + self.l2*cos(theta2)
        return (x1, y1), (x2, y2)
    
    def step(self, t, dt):
        self.solution += self.RK4(t, self.solution, dt)