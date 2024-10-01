import colorsys
import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__))) # to fix relative import errors

import numpy
import pygame
import sys

from numpy import sin, cos
from scipy.constants import g, pi
from pygame import gfxdraw
from multipendulum import Pendulum

# constants
# 0.7295790747075395, 0.6918422260157702
# pi, pi, pi - 0.001
start_thetas = [0.7295790747075395, 0.6918422260157702, 0]
pendulum = Pendulum(start_thetas, [1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
masses, lengths = pendulum.masses, pendulum.lengths

N = pendulum.n
SCALE = 400/N
FPS = 120
resolution = (850, 850) # FIXED, DO NOT CHANGE
pivot_position = tuple(d/2 for d in resolution)
zoom = False

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
TRAIL_COLOR = (230, 230, 230)
BOB_COLOR = (0, 0, 200)

# functions

def graphics(bobs, previous_point, trace_screen, trail_color = TRAIL_COLOR):
    scale = 15/N
    X, Y = [int(bob[0]) for bob in bobs], [int(bob[1]) for bob in bobs]
    
    if previous_point:
        pygame.draw.aaline(trace_screen, trail_color, (previous_point[0], previous_point[1]), bobs[N - 1], 3)
        
    screen.fill(WHITE)
    screen.blit(trace_screen, (0, 0))
    
    # Draw Pivot
    gfxdraw.aacircle(screen, int(pivot_position[0]), int(pivot_position[1]), int(scale/2), BLACK) # pivot outline
    gfxdraw.filled_circle(screen, int(pivot_position[0]), int(pivot_position[1]), int(scale/2), BLACK) # pivot outline
    
    # Pivot Rod
    pygame.draw.aaline(screen, BLACK, pivot_position, bobs[0], 3)
    
    # Draw Rods
    for n in range(0, N - 1):
        pygame.draw.aaline(screen, BLACK, bobs[n], bobs[n + 1], 3)
    
    # Draw Bobs
    for n in range(0, N):
        gfxdraw.aacircle(screen, X[n], Y[n], int(masses[n]*scale), BLACK) # bob1 outline
        gfxdraw.filled_circle(screen, X[n], Y[n], int(masses[n]*scale), BOB_COLOR) # bob1 fill
        
    return bobs[N - 1]

def fps_counter(screen):
    fps_t = font.render(str(int(clock.get_fps())) + " FPS", 1, pygame.Color("RED"))
    screen.blit(fps_t,(0, 0))
    
def debug_info(screen, t, tc):
    debug_top = fontDebug.render(f'{N} pendulums - time elapsed: {str(round(t, 1))} seconds, dt: {round(1/FPS, 4)}, masses: {masses}, lengths: {lengths}', 100, BLACK)
    debug_lower = fontDebug.render(f'Number of trace lines: {tc}', 100, BLACK)
    screen.blit(debug_top, (100, 0))
    screen.blit(debug_lower, (100, 20))

# pygame inits
pygame.init()
screen = pygame.display.set_mode(resolution)
pygame.display.set_caption(f'{N} Pendulum')
screen.fill(WHITE)
pygame.display.update()
clock = pygame.time.Clock()

font = pygame.font.SysFont("songti", 20)
fontDebug = pygame.font.SysFont("songti", 20 - int(N))

def animate():
    # Vars
    dt = 1/FPS
    t = 0
    trace_count = -1

    trace_screen = screen.copy()
    previous_point = []

    trace_screen.fill(WHITE)

    running = True

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        bobs = pendulum.position(scale = (-SCALE, SCALE), shift = pivot_position)
        
        previous_point = graphics(bobs, previous_point, trace_screen)
        trace_count += 1

        t += dt
        
        pendulum.step(t, dt)
        
        fps_counter(screen)
        debug_info(screen, t, trace_count)
        
        clock.tick(FPS)
        pygame.display.update()

animate()