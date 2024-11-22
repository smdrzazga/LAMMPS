from common.SCREEN import *
from common.CONTAINERS import *
import numpy as np
from scipy.optimize import curve_fit

class NTBPhase:
    def __init__(self, periods) -> None:
        self.periods = periods
        self.C_left = 0 + 0j
        self.C_right = 0 + 0j
        self.C = 0 + 0j

    def clear(self):
        self.C_left = 0 + 0j
        self.C_right = 0 + 0j
        self.C = 0 + 0j

    def update(self, molecule: Banana, box: SimulationBox) -> None:
        center = molecule.center_of_mass()
        polarization = molecule.polarization()[1]

        self.C += polarization * np.exp(2j*self.periods*np.pi * center[2] / box.get_side_length(2))
        if center[0] < box.get_side_length(0)//2:
            self.C_left += polarization * np.exp(2j*self.periods*np.pi * center[2] / box.get_side_length(2))
        else:
            self.C_right += polarization * np.exp(2j*self.periods*np.pi * center[2] / box.get_side_length(2))

    def compute_com_drift(self, box: SimulationBox) -> float:
        drift = box.get_side_length(2) / self.periods * np.angle(self.C) / (2*np.pi)
        return drift


class SmecticParameter:
    def __init__(self, SMECTIC_PERIODS: int) -> None:
        self.parameter = 0 + 0j
        self.count = 0
        self.SMECTIC_PERIODS = SMECTIC_PERIODS

    def __repr__(self) -> str:
       return f"{np.abs(self.parameter):.3f} | {self.count}"

    def calculate_smectic_param(self, screen: np.array) -> float:
        y_max, x_max = screen.shape
        for x in range(x_max):
            for y in range(y_max):
                self.parameter += np.exp(self.SMECTIC_PERIODS * 2*np.pi*1j * y / y_max) * screen[y, x]
                self.count += screen[y, x]

        if self.count != 0:
            return np.absolute(self.parameter) / self.count
        return 0.0    
    

class CorrelationLength:
    def _decay_func(self, x, A, L) -> float:
        return A*np.exp(-1*x / L)

    def _make_symmetrized_range(self, data: np.array) -> np.array:
        max_coord = data[0] + data[-1]
        sym_data = data
        for i in range(len(sym_data)):
            if sym_data[i] > max_coord / 2:
                sym_data[i] -= 2*(sym_data[i] - max_coord/2)
        return sym_data

    def fit_decay(self, data) -> list[float, float, float, float]:
        x_data = data[2:-2, 0]
        x_data = self._make_symmetrized_range(x_data)
        y_data = data[2:-2, 1]

        popt, pvar = curve_fit(self._decay_func, x_data, y_data, p0=(0.1, 0.5))
        result = [popt[0], np.sqrt(pvar[0][0]), popt[1], np.sqrt(pvar[1][1])]
        return result


class EllipsisSemiaxes:
    def calculate_semiaxes(self, slice):
            x = slice[:, :, 0].flatten()
            x -= np.average(x)

            y = slice[:, :, 1].flatten()
            y -= np.average(y)

            A = np.stack([x**2, y**2]).T
            b = np.ones_like(x)
            w = np.linalg.lstsq(A, b, rcond=-1)[0].squeeze()
                        
            return np.sqrt(1/w)