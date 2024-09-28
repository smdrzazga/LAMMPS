import numpy as np

class HeatmapReader:
    def read_matrix(file: str, n_slices: float, n_periods: int) -> list:
        list_of_params = np.zeros(n_slices, dtype=np.complex128)
        list_of_counts = np.zeros(n_slices)

        with open(file, "r") as f:
            for pixel in f:
                # read colour and coordinates of pixels
                z, y, x, N = pixel.split()
                z, y, x, N = float(z), float(y), float(x), int(N)

                # calculate slice to which pixel belongs and accumulate its contribution to smectic parameter
                index = int(x * n_slices)
                list_of_params[index] += np.exp(n_periods* 2*np.pi*1j * y) * N
                list_of_counts[index] += N

        for i in range(n_slices):
            if list_of_counts[i] != 0:
                list_of_params[i] /= list_of_counts[i]

        return list_of_params

