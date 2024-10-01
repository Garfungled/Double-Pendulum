import numpy

from matplotlib import cm
from matplotlib import pyplot
from matplotlib.ticker import LinearLocator
from scipy.constants import pi

def read_heatmap_from_csv(file):
    return numpy.loadtxt(file, delimiter=',')

def plot_3dheatmap(heatmap, theta1_range, theta2_range, theta3_range):
    fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
    fig.set_size_inches(10, 10)

    # Make data.
    theta1_range, theta2_range = numpy.meshgrid(theta1_range, theta2_range)
    Z = heatmap

    # Plot the surface.
    surf = ax.plot_surface(theta1_range, theta2_range, Z, cmap="afmhot", linewidth=0, antialiased=False, norm="symlog")
    
    # Configure Heatmap axis
    pyplot.colorbar(surf, ticks=[0.5, 1, 1.5, 2, 3, 4], label='Lyapunov Exponent').ax.set_yticklabels(['0.5', '1', '1.5', '2', '3', '4'])
    
    # Configure Axes
    pyplot.xlabel("Theta 2 (radians)")
    pyplot.xticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.50*pi, 1.75*pi, 2*pi], ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'])
        
    pyplot.ylabel("Theta 1 (radians)")
    pyplot.yticks([0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.50*pi, 1.75*pi, 2*pi], ['0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'])

    ax.set_zlim(0, numpy.max(Z))
    ax.zaxis.set_major_locator(LinearLocator(10))
    
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.

    pyplot.show()
    

heatmap_file = "data/500x500/data_0.csv"
heatmap = read_heatmap_from_csv(heatmap_file)
dimension = len(heatmap[0])

theta1_range, theta2_range, theta3_range = numpy.linspace(0, 2 * pi, dimension), numpy.linspace(0, 2 * pi, dimension), numpy.linspace(0, 2 * pi, dimension)

plot_3dheatmap(heatmap, theta1_range, theta2_range, theta3_range)