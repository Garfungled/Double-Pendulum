import sys
from time import sleep
import numpy

from matplotlib import pyplot
from scipy.constants import pi

# file looks like "num1,num2,num3,num4,num5,..."
def read_heatmap_from_csv(file):
    return numpy.loadtxt(file, delimiter=',')

def on_window_close(event):
    sys.exit()

def plot_heatmap(heatmaps, theta1_range, theta2_range, animate=False):
    dimension = len(heatmaps[0])
    rTp = lambda n : int(n / 2/pi * dimension)
    
    figure = pyplot.figure(figsize = (10, 10))
    extent = [theta2_range.min(), theta2_range.max(), theta1_range.min(), theta1_range.max()]
    axesimage = pyplot.imshow(heatmaps[0], origin='lower', cmap='hot', extent=extent, aspect='equal', norm="symlog")
    
    pyplot.colorbar(axesimage, ticks=[0, 1, 2, 3, 4, 5, 6, 7], label='Lyapunov Exponent').ax.set_yticklabels(['0', '1', '2', '3', '4', '5', '6', '7'])
    pyplot.xlabel("Theta 2 (radians)")
    pyplot.xticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.50*pi, 1.75*pi, 2*pi], ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'])    
    pyplot.ylabel("Theta 1 (radians)")
    pyplot.yticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.50*pi, 1.75*pi, 2*pi], ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'])
    pyplot.tight_layout()
    
    def coords(x, y):
        if x >= theta2_range.min() and x <= theta2_range.max() and y >= theta1_range.min() and y <= theta1_range.max():
            xi, yi = rTp(x), rTp(y)
            return f'(Theta 1 = {theta1_range[yi]}, Theta 2 = {theta2_range[xi]})(x = {xi}, y = {yi})'
    
    def on_key(event):
        if event.key == 'a':
            x, y = event.xdata, event.ydata
            if x is not None and y is not None and theta2_range.min() <= x <= theta2_range.max() and theta1_range.min() <= y <= theta1_range.max():
                # Do something with the data
                print(f"{theta1_range[rTp(x)]},{theta2_range[rTp(y)]}")
    
    pyplot.gcf().canvas.mpl_connect('key_press_event', on_key)
    pyplot.gca().format_coord = coords
    
    if animate:
        while True:
            figure.canvas.mpl_connect('close_event', on_window_close)
            
            i = 0
            for heatmap in heatmaps:
                axesimage.set_data(heatmap)
                pyplot.draw()
                pyplot.title(f"{i/(step) * 2}π")
                pyplot.pause(1/60)
                i += 1
    else:
        pyplot.show()


dimension = 1000
step = 64

# heatmap_file = f"data/data_{dimension}x{dimension}.csv"
# heatmap = read_heatmap_from_csv(heatmap_file)

theta1_range, theta2_range = numpy.linspace(0, 2 * pi, dimension), numpy.linspace(0, 2 * pi, dimension)
heatmap = [read_heatmap_from_csv(f"data/{dimension}x{dimension}/data_{i}.csv") for i in range(0, step)]

plot_heatmap(heatmap, theta1_range, theta2_range, animate=True)