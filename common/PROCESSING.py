from main import GlobalParameters
from multiprocessing.pool import Pool
from concurrent.futures import ProcessPoolExecutor

from IO import *
from SCREEN import *
from CONTAINERS import *
from ORDER_PARAMETERS import *
from TIMERS import *

class BatchAnalyzer:
    def __init__(self, parameters: GlobalParameters) -> None:
        self.proc_params = parameters.proc_params
        self.file_params = parameters.file_params
        self.snap_params = parameters.snap_params
        self.reader = LAMMPSReader(parameters.file_params)
        self.binner = AtomBinner(self.snap_params['size'], view_plane=self.snap_params['plane'])
        self.phase_tracker = NTBPhase(self.snap_params['DIRECTOR_PERIODS'])
        self.screen = None
        self.screenshot = None

    def setup(self):
        boundaries = self.reader.read_boundaries()
        N_ATOMS = self.reader.read_number_of_atoms()
        self.box = SimulationBox(boundaries, N_ATOMS)
        self.molecule = Molecule(1, self.snap_params['ATOMS_IN_MOLECULE'])
        self.atom = Atom()

    def analyze_batch(self, ID) -> Screen:
        self.reader.open()
        self.setup()

        line = self.reader.get_line_split()
        while line:
            try:
                self.read_atom(line)
            except:
                continue
                     
            self.molecule.add(self.atom)

            if not self.molecule.is_full():
                continue

            if self.snap_params['ANALYSE_WALL'] and not self.molecule.is_at_wall():
                continue

            self.add_molecule_to_pixel()
            self.molecule.clear()

            if self.atom.is_last(self.box.get_num_atoms()):
                if self.snap_params['ELIMINATE_GOLDSTONE']:
                    self._eliminate_goldstone_mods_both_sides()
                self.finalize_screenshot_analysis()
                
            line = self.reader.get_line_split()

        print(f"Task is done!: {ID}")
        return self.screen

    def analyze_batch(self, ID) -> NotImplementedError:
        raise NotImplementedError("Function 'analyze_batch' is virtual in this scope.")

    def add_molecule_to_pixel(self) -> NotImplementedError:
        raise NotImplementedError("Function 'add_molecule_to_pixel' is virtual in this scope.")
    
    def _eliminate_goldstone_mods_both_sides(self) -> NotImplementedError:
        raise NotImplementedError("Function '_eliminate_goldstone_mods_both_sides' is virtual in this scope.")

    def clear_screenshots(self) -> NotImplementedError:
        raise NotImplementedError("Function 'clear_screenshots' is virtual in this scope.")
    
    def read_atom(self, split_line):
        self.atom = Atom(split_line[-3:], id=split_line[0])

    def finalize_screenshot_analysis(self):
        self.screen.append_screenshot(self.screenshot)
        self.update_box()
        self.clear_screenshots()
        self.phase_tracker.clear()
    
    def update_box(self):
        boundaries = self.reader.read_boundaries()
        self.box.update_boundaries(boundaries)


class DirectorFullAnalyzer(BatchAnalyzer):
    def __init__(self, parameters: GlobalParameters) -> None:
        super().__init__(parameters)
        size = self.snap_params['size']
        
        self.screen = Screen(size, DirectorPixel)
        self.screenshot = Screenshot(size, DirectorPixel)
        self.screenshotCenter = Screenshot(size, CenterPixel)
    
    def setup(self):
        boundaries = self.reader.read_boundaries()
        N_ATOMS = self.reader.read_number_of_atoms()
        self.box = SimulationBox(boundaries, N_ATOMS)
        self.molecule = Banana(1, self.snap_params['ATOMS_IN_MOLECULE'])
        self.atom = Atom()

    def add_molecule_to_pixel(self):
        self.phase_tracker.update(self.molecule, self.box)

        director = Atom(self.molecule.director())
        center = Atom(self.molecule.center_of_mass())
        center.wrap(self.box.get_all_side_lengths())
        pixel_position = self.binner.determine_pixel(center)

        self.screenshot.assign(director, pixel_position)
        self.screenshotCenter.assign(center, pixel_position)

    def clear_screenshots(self):
        self.screenshot = Screenshot(self.snap_params['size'], DirectorPixel)
        self.screenshotCenter = Screenshot(self.snap_params['size'], CenterPixel)

    def _eliminate_goldstone_mods_both_sides(self, tracker: NTBPhase, box: SimulationBox) -> None:
        drift = tracker.compute_com_drift(box)
        pix_to_scroll_both = self.screenshotCenter.pixels_to_scroll(self.snap_params['size'][1], box, drift)
        self.screenshot.scroll_both_sides(pix_to_scroll_both)

    def get_results_in_parallel(self) -> None:
        with Pool(self.proc_params['NP']) as executor:
            IDs = [ID for ID in range(self.proc_params['NP'])]
            for result in executor.starmap(self.analyze_batch, IDs):
                self.screen.append_screenshot(result)


