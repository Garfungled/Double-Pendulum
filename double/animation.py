import colorsys
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__))) # to fix relative import errors

import numpy
import pygame
import sys

from numpy import sin, cos
from scipy.constants import g, pi
from pygame import gfxdraw
from pendulum import Pendulum

# constants

m1, m2 = 1.0, 1.0
l1, l2 = 1.0, 1.0

FPS = 120
resolution = (850, 850)
pivot_position = tuple(d/2 for d in resolution)
zoom = False

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
TRAIL_COLOR = (230, 230, 230)
BOB1_COLOR = (0, 0, 230)
BOB2_COLOR = (0, 0, 230)

# functions

def F(t, y):
    theta1, theta2, omega1, omega2 = y
    alpha1 = ((m2*g*sin(theta2)*cos(theta1 - theta2)) - (m2*sin(theta1 - theta2))*(l1*(omega1**2)*cos(theta1 - theta2) + l2*(omega2**2)) - g*sin(theta1)*(m1 + m2))/(l1*(m1 + m2*(sin(theta1 - theta2)**2)))
    alpha2 = ((m1 + m2)*((l1*(omega1**2)*sin(theta1 - theta2)) - (g*sin(theta2)) + (g*sin(theta1)*cos(theta1 - theta2))) + (m2*l2*(omega2**2)*sin(theta1 - theta2)*cos(theta1 - theta2)))/(l2*(m1 + m2*(sin(theta1 - theta2)**2)))

    return numpy.array([omega1, omega2, alpha1, alpha2])

def RK4(t, y, dt):
    k1 = F(t, y)
    k2 = F(t + 0.5*dt, y + 0.5*k1*dt)
    k3 = F(t + 0.5*dt, y + 0.5*k2*dt)
    k4 = F(t + dt, y + dt*k3)
    
    return (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

def update(theta1, theta2):
    scale = 200 * (l1 + l2)/2

    x1 = l1*scale*sin(theta1) + pivot_position[0]
    y1 = l1*scale*cos(theta1) + pivot_position[1]
    x2 = x1 + l2*scale*sin(theta2)
    y2 = y1 + l2*scale*cos(theta2)

    return (x1, y1), (x2, y2)

def adjust(n):
    pass

def graphics(bob1, bob2, previous_point, trace_screen, screen, trail_color = TRAIL_COLOR):
    scale = 7
    x1, y1 = int(bob1[0]), int(bob1[1])
    x2, y2 = int(bob2[0]), int(bob2[1])

    if previous_point:
        pygame.draw.aaline(trace_screen, trail_color, (previous_point[0], previous_point[1]), (x2, y2), 3)

    screen.fill(WHITE)
    screen.blit(trace_screen, (0, 0))

    pygame.draw.aaline(screen, BLACK, pivot_position, bob1, 3) # rod1
    pygame.draw.aaline(screen, BLACK, bob1, bob2, 3) # rod1

    gfxdraw.aacircle(screen, int(pivot_position[0]), int(pivot_position[1]), 5, BLACK) # pivot outline
    gfxdraw.filled_circle(screen, int(pivot_position[0]), int(pivot_position[1]), 5, BLACK) # pivot outline

    gfxdraw.aacircle(screen, x1, y1, int(m1*scale), BLACK) # bob1 outline
    gfxdraw.filled_circle(screen, x1, y1, int(m1*scale), BOB1_COLOR) # bob1 fill

    gfxdraw.aacircle(screen, x2, y2, int(m2*scale), BLACK) # bob2 outline
    gfxdraw.filled_circle(screen, x2, y2, int(m2*scale), BOB2_COLOR) # bob2 fill

    return bob2

def fps_counter(screen, clock):
    fps_t = font.render(str(int(clock.get_fps())) + " FPS", 1, pygame.Color("RED"))
    screen.blit(fps_t,(0,0))
    
def debug_info(screen, t, tc):
    debug_top = fontDebug.render(f'2 pendulums - time elapsed: {str(round(t, 1))} seconds, dt: {round(1/FPS, 4)}, masses: {(m1, m2)}, lengths: {(l1, l2)}', 100, BLACK)
    debug_lower = fontDebug.render(f'Number of trace lines: {tc}', 100, BLACK)
    screen.blit(debug_top, (100, 0))
    screen.blit(debug_lower, (100, 20))

# pygame initialization

pygame.init()
font = pygame.font.SysFont("Arial", 20)
fontDebug = pygame.font.SysFont("songti", 18)

def animate(theta1, theta2):
    pygame.init()
    screen = pygame.display.set_mode(resolution)
    pygame.display.set_caption("Single Pendulum")
    screen.fill(WHITE)
    pygame.display.update()
    clock = pygame.time.Clock()
    
    dt = 1/FPS
    t = 0
    trace_count = -1

    trace_screen = screen.copy()
    previous_point = []

    solution = numpy.array([theta1, theta2, 0, 0]) # theta1, theta2, omega1, omega2

    trace_screen.fill(WHITE)

    running = True

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        bob1, bob2 = update(solution[0], solution[1])
        
        previous_point = graphics(bob1, bob2, previous_point, trace_screen, screen)
        trace_count += 1

        t += dt

        solution += RK4(t, solution, dt)
        
        fps_counter(screen, clock)
        debug_info(screen, t, trace_count)
        
        clock.tick(FPS)
        pygame.display.update()
        
animate(3.220937704355627,6.105640936159007)