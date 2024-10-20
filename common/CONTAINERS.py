from common.MOLECULES import *
from config import *
from mmap import ALLOCATIONGRANULARITY
import os

class SimulationData:
    def __init__(self, location) -> None:
        self.location = location
        self.size = self.get_file_size()

    def get_location(self):
        return self.location
    
    def get_file_size(self):
        file_stats = os.stat(self.location)
        file_size = file_stats.st_size
        return file_size

class ChunkData(SimulationData):
    def __init__(self, proc_params: ProcessingParameters) -> None:
        super().__init__(proc_params.params["INPUT_FILE"])
        self.proc_params = proc_params.params

    def get_offset(self, ID) -> int:
        start_percent, stop_percent = self.proc_params["ANALYZE_RANGE"]
        offset = int(self.get_file_size()*start_percent + self.get_chunk_size()*ID)
        return offset
    
    def get_chunk_size(self) -> int:
        file_percent = self.proc_params["ANALYZE_RANGE"][1] - self.proc_params["ANALYZE_RANGE"][0]
        chunk_size = self.get_file_size() * file_percent // self.proc_params["NP"]  
        return self._round_up(chunk_size, ALLOCATIONGRANULARITY)
    
    def _round_up(self, value, precision) -> int:
        value += precision
        value -= value % precision
        return int(value)
    

class SimulationBox:
    def __init__(self, boundaries: list, n_atoms) -> None:
        self.update_boundaries(boundaries)

        self.volume = np.prod(self.size)
        self.atoms = int(n_atoms)

    def get_all_side_lengths(self) -> list:
        return self.size

    def get_side_length(self, i) -> float:
        return self.size[i]
    
    def get_num_atoms(self):
        return self.atoms
    
    def get_volume(self):
        return self.volume
    
    def update_boundaries(self, new_boundaries: list):
        [x_min, x_max, y_min, y_max, z_min, z_max] = new_boundaries
        self.min = np.array([x_min, y_min, z_min], dtype=np.float32)
        self.max = np.array([x_max, y_max, z_max], dtype=np.float32)
        self.size = np.array([x_max - x_min, y_max - y_min, z_max - z_min], dtype=np.float32)

