from CONTAINERS import *
import numpy as np

class NTBPhaseTracker:
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

    def compute_com_drift(self, box: SimulationBox):
        drift = box.get_side_length(2) / self.periods * np.angle(self.C) / (2*np.pi)
        return drift
