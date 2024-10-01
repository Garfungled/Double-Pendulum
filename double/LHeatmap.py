import numpy

from matplotlib import pyplot
from scipy.constants import pi

# file looks like "num1,num2,num3,num4,num5,..."
def read_heatmap_from_csv(file):
    return numpy.loadtxt(file, delimiter=',')

def plot_heatmap(heatmap, theta1_range, theta2_range):
    dimension = len(heatmap)
    pyplot.figure(figsize = (10, 10))
    
    extent = [theta2_range.min(), theta2_range.max(), theta1_range.min(), theta1_range.max()]
    
    axesimage = pyplot.imshow(heatmap, origin='lower', cmap='hot', extent=extent, aspect='equal', norm="symlog")
    pyplot.colorbar(axesimage, ticks=[0, 1, 2, 3, 4, 5], label='Lyapunov Exponent').ax.set_yticklabels(['0', '1', '2', '3', '4', '5'])
    
    pyplot.xlabel("Theta 2 (radians)")
    pyplot.xticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.50*pi, 1.75*pi, 2*pi], ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'])    
    pyplot.ylabel("Theta 1 (radians)")
    pyplot.yticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.50*pi, 1.75*pi, 2*pi], ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'])
    
    def coords(x, y):
        if x >= theta2_range.min() and x <= theta2_range.max() and y >= theta1_range.min() and y <= theta1_range.max():
            xi = int(x / 2/pi * dimension)
            yi = int(y / 2/pi * dimension)
            return f'(Theta 1 = {theta1_range[yi]}, Theta 2 = {theta2_range[xi]})(x = {xi}, y = {yi})'
        else:
            return ''

    pyplot.gca().format_coord = coords
    
    theta1 = None
    theta2 = None
    def on_key(event):
        global theta1, theta2
        if event.key == 'a':
            x, y = event.xdata, event.ydata
            if x is not None and y is not None:
                if theta2_range.min() <= x <= theta2_range.max() and theta1_range.min() <= y <= theta1_range.max():
                    xi = int(x / (2 * pi) * dimension)
                    yi = int(y / (2 * pi) * dimension)
                    theta1 = theta1_range[yi]
                    theta2 = theta2_range[xi]
                    # Do something with the data
                    print(f"{theta1},{theta2}")
    
    pyplot.gcf().canvas.mpl_connect('key_press_event', on_key)
    
    pyplot.tight_layout()
    pyplot.show()
    
    
heatmap_file = "fastData/data/data_2000x2000.csv"
heatmap = read_heatmap_from_csv(heatmap_file)

dimension = len(heatmap[0])

theta1_range, theta2_range = numpy.linspace(0, 2 * pi, dimension), numpy.linspace(0, 2 * pi, dimension)

plot_heatmap(heatmap, theta1_range, theta2_range)