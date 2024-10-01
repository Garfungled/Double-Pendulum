import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__))) # to fix relative import errors

import csv
import unicodedata
import lyapunov as lyapunov

from scipy.constants import pi

dimension = 300

# from 0 to 2pi for each
with open(f"data/data_5.csv", mode = "w", newline = '') as file:
    writer = csv.writer(file)
    
    for row in range(0, dimension):
        for col in range(0, dimension):
            file.write(str(lyapunov.lyapunov((2 * pi)/dimension * row, (2 * pi)/dimension * col, 0.001)))
            file.write('' if col == dimension - 1 else ',')
        file.write("\n")