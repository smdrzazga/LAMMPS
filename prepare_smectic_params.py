import numpy as np
import os
from common.SCREEN import *
from common.ORDER_PARAMETERS import *
from scipy.optimize import curve_fit


class CentersMatrix:
    def __init__(self, size: tuple[int, int] = (150, 150)) -> None:
        self.size = size

    def read_matrix_from_file(self, file):
        with open(file, 'r') as f:
            data = np.array([line.split()[-1] for line in f], dtype=np.float16)
        self.matrix = data.reshape(self.size)

    def update(self, new_matrix: np.array) -> None:
        self.matrix = new_matrix

    def get_slice(self, x_lim = (None, None), y_lim = (None, None)) -> np.array:
        return self.matrix[y_lim[0]:y_lim[1], x_lim[0]:x_lim[1]]



if __name__ == "__main__":
    file = "C:/Users/" + os.getlogin() + "/Desktop/LAMMPS_matrices/centers_matrices/centers_screen_bulk_6k_0.32.txt"

    N = 150
    N_SLICES = 30
    SMECTIC_PERIODS = 4

    matrix = CentersMatrix(size=(N, N))
    smec = SmecticParameter(periods = SMECTIC_PERIODS)

    matrix.read_matrix_from_file(file)
    slice = matrix.get_slice(x_lim = (2,8))

    print(smec.calculate_smectic_param(slice))
