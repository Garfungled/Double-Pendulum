import numpy
from numpy import sin, cos
from scipy.constants import g, pi

class Pendulum:
    def __init__(self, thetas: list, lengths: tuple, masses: tuple):
        self.n = len(thetas)
        self.lengths = lengths
        self.masses = masses
        
        self.solution = numpy.array([[-theta for theta in thetas], [0] * len(thetas)])
        
    
    def __F(self, t, y):
        n = self.n
        thetas, omegas = y[0], y[1]
        sigma = lambda i, j: i if j <= i else j
        
        A = []
        for i in range(0, n):
            A.append([])
            for j in range(0, n):
                A[i].append(self.lengths[j] * cos(thetas[i] - thetas[j]) * sum([self.masses[k] for k in range(sigma(i, j), n)]))
                
        B = []
        for i in range(0, n):
            B.append(-g * sin(thetas[i]) * sum([self.masses[j] for j in range(i, n)]) - sum([self.lengths[j] * (omegas[j]**2) * sin(thetas[i] - thetas[j]) * sum([self.masses[k] for k in range(sigma(i, j), n)]) for j in range(0, n)]))
        
        X = list(numpy.linalg.inv(A).dot(B)) # Alphas
        
        return numpy.array([omegas, X])
        
    
    def RK4(self, t, y, dt):
        k1 = self.__F(t, y)
        k2 = self.__F(t + 0.5*dt, y + 0.5*k1*dt)
        k3 = self.__F(t + 0.5*dt, y + 0.5*k2*dt)
        k4 = self.__F(t + dt, y + dt*k3)
        
        return (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    def position(self, scale = (1, 1), shift = (0, 0)):
        thetas = self.solution[0]
        X = [sum([self.lengths[k] * sin(thetas[k]) for k in range(0, n)]) * scale[0] + shift[0] for n in range(1, self.n + 1)]
        Y = [-sum([self.lengths[k] * cos(thetas[k]) for k in range(0, n)]) * scale[0] + shift[0] for n in range(1, self.n + 1)]
        
        return [(X[n], Y[n]) for n in range(0, self.n)]
        
    def step(self, t, dt):
        self.solution += self.RK4(t, self.solution, dt)