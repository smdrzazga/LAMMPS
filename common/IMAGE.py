from common.IO import Writer
from common.SCREEN import Screen
from common.ORDER_PARAMETERS import *
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


class ScatterPrinter:
    def __init__(self) -> None:
        self.scatter = None

    def create_scatter(self, screen: Screen, N_SLICES) -> None:
        slices = self._create_slices(screen, N_SLICES)
        slice_directors = self.calculate_directors(slices)
        return slice_directors[:2, :]

    def _create_slices(self, screen: Screen, N_SLICES) -> list[Slice]:
        slice_width = screen.get_size_in_pixels()[0] // N_SLICES
        sizes = [[[slice_width*(i-1), slice_width*i], [None, None]] for i in range(0, N_SLICES)]
        slices = [Slice(screen, sizes[i]) for i in range(0, N_SLICES)]
        return slices
    
    def calculate_directors(self, slices: list[Slice]) -> list[list]:
        slice_directors = np.zeros((len(slices), 3))
        for i, slice in enumerate(slices):
            slice_directors[i] = slice.slice_director()
        return slice_directors

    def save_scatter(self, data, target_location) -> None:
        writer = Writer(target_location)
        writer.clear()
        writer.open()
        for i, point in enumerate(data):
            line = f"{i} "
            for coord in point:
                line += f"{coord} "
            writer.write_line(point + '\n')
        writer.close()