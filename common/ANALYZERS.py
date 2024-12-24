from common.IO import *
from common.SCREEN import *
from common.CONTAINERS import *
from common.ORDER_PARAMETERS import *
from common.TIMERS import *

from config import GlobalParameters
from multiprocessing.pool import Pool
import glob


class BatchAnalyzer:
    def __init__(self, parameters: GlobalParameters) -> None:
        self.proc_params = parameters.proc_params
        self.file_params = parameters.file_params
        self.snap_params = parameters.snap_params
        self.reader = LAMMPSReader(self.proc_params)
        self.binner = AtomBinner(self.snap_params['SIZE'], view_plane=self.snap_params['PLANE'])
        self.phase_tracker = NTBPhase(self.snap_params['DIRECTOR_PERIODS'])
        self.screen: Screen
        self.screenshot: Screenshot

    def setup(self):
        boundaries = self.reader.read_boundaries()
        N_ATOMS = self.reader.read_number_of_atoms()
        self.box = SimulationBox(boundaries, N_ATOMS)
        self.molecule = Molecule(1, self.snap_params['ATOMS_IN_MOLECULE'])
        self.atom = Atom()

    def analyze_batch(self, ID) -> NotImplementedError:
        raise NotImplementedError("Function is virtual in this scope")
    
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

    def clear_screenshots(self) -> None:
        raise NotImplementedError("Function is virtual in this scope.")


class DirectorFullAnalyzer(BatchAnalyzer):
    def __init__(self, parameters: GlobalParameters) -> None:
        super().__init__(parameters)
        size = self.snap_params['SIZE']
        
        self.screen = Screen(size, DirectorPixel)
        self.screenshot = Screenshot(size, DirectorPixel)
        self.screenshotCenter = Screenshot(size, CenterPixel)
    
    def analyze_batch(self, ID) -> Screen:
        self.reader.open(ID)
        self.setup()

        line = self.reader.get_line_split()
        while line:
            try:
                self.read_atom(line)
            except:
                line = self.reader.get_line_split()
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

    def setup(self):
        N_ATOMS = self.reader.read_number_of_atoms()
        boundaries = self.reader.read_boundaries()
        self.box = SimulationBox(boundaries, N_ATOMS)
        self.molecule = Banana(1, self.snap_params['ATOMS_IN_MOLECULE'])
        self.atom = Atom()

    def add_molecule_to_pixel(self):
        self.phase_tracker.update(self.molecule, self.box)

        director = Atom(self.molecule.director())
        center = Atom(self.molecule.center_of_mass())
        center.wrap(self.box.get_all_side_lengths())
        pixel_position = self.binner.determine_pixel(center, self.box)

        self.screenshot.assign(director, pixel_position)
        self.screenshotCenter.assign(center, pixel_position)

    def clear_screenshots(self):
        self.screenshot = Screenshot(self.snap_params['SIZE'], DirectorPixel)
        self.screenshotCenter = Screenshot(self.snap_params['SIZE'], CenterPixel)

    def _eliminate_goldstone_mods_both_sides(self, tracker: NTBPhase, box: SimulationBox) -> None:
        drift = tracker.compute_com_drift(box)
        pix_to_scroll_both = self.screenshotCenter.pixels_to_scroll(self.snap_params['SIZE'][1], box, drift)
        self.screenshot.scroll_both_sides(pix_to_scroll_both)

    def get_results_in_parallel(self) -> Screen:
        with Pool(self.proc_params['NP']) as executor:
            IDs = zip([ID for ID in range(self.proc_params['NP'])])
            for result in executor.starmap(self.analyze_batch, IDs):
                self.screen.append_screenshot(result)

        return self.screen
    
    def get_results_in_series(self) -> None:
        IDs = [ID for ID in range(self.proc_params['NP'])]
        for ID in IDs:
            result = self.analyze_batch(ID)
            self.screen.append_screenshot(result)

        return self.screen


class ParameterAnalyzer:
    def __init__(self):
        self.set_source_paths()

    def set_source_paths(self) -> str:
        self.matrices_dir    = "C:\\Users\\" + os.getlogin() + "\\Desktop\\LAMMPS_matrices"
        self.smectics_dir    = self.matrices_dir + "\\smectic_params"
        self.centers_dir     = self.matrices_dir + "\\centers_matrices"
        self.directors_dir   = self.matrices_dir + "\\directors_matrices"
        self.correlation_dir = self.matrices_dir + "\\correlation_lengths"
        self.ellipsis_dir    = self.matrices_dir + "\\ellipsis_semiaxes"

    def get_all_files_in_dir_to_analyze(self, dir: str, filename_prefix: str) -> list[str]:
        files_to_analyze = glob.glob(dir + '\\' + filename_prefix)
        return files_to_analyze
  
    def print_results_to_file(self, results: list[tuple[str, ...]], file: str) -> None:
        with open(file, 'w+') as t:
            for line in results:
                print(' '.join(line), file=t)

    def get_density_from_path(self, filepath) -> str:
        return filepath.split('_')[-1][:-4]   


class SmecticAnalyzer(ParameterAnalyzer):
    def __init__(self, N_SLICES: int, SMECTIC_PERIODS: int, N_PIX: int = 150) -> None:
        self.N_SLICES = N_SLICES
        self.SMECTIC_PERIODS = SMECTIC_PERIODS
        self.N = N_PIX
        self.set_source_paths()


    def print_smectic_params_for_all_files(self, centers_file_prefix: str) -> None:
        files_to_analyze = self.get_all_files_in_dir_to_analyze(self.centers_dir, centers_file_prefix)
        print("Files to analyze: ", files_to_analyze)
        for file in files_to_analyze:
            self.print_smectic_params_for_single_file(file)


    def print_smectic_params_for_single_file(self, file: str) -> None:
        matrix = self.read_centers_data(file)

        parameters = self.calculate_smectic_params_from_matrix(matrix)
        results = self.prepare_results_to_print(parameters)
        target = self.target_location(file)
        self.print_results_to_file(results, target)


    def read_centers_data(self, file) -> DirectorsMatrix:
        matrix = CentersMatrix((self.N, self.N))
        matrix.read_matrix_from_file(file)
        return matrix


    def calculate_smectic_params_from_matrix(self, matrix: CentersMatrix) -> list[float]:
        screen_width = matrix.size[-1]
        step_size = screen_width // self.N_SLICES
        parameters = []

        for i in range(0, screen_width, step_size):
            slice = matrix.get_slice(x_lim=(i, i + step_size))
            parameter = self.calculate_smectic_param_for_slice(slice)
            parameters.append(parameter)
        return parameters


    def calculate_smectic_param_for_slice(self, slice) -> float:
        return SmecticParameter(SMECTIC_PERIODS = self.SMECTIC_PERIODS).calculate_smectic_param(slice)
        

    def prepare_results_to_print(self, parameters: list[float]) -> list[tuple[str, str]]:
        slice_positions = [str((i+1/2)/self.N_SLICES) for i in range(self.N_SLICES)]
        str_parameters = [str(param) for param in parameters]
        results = zip(slice_positions, str_parameters)
        return results


    def target_location(self, source_path) -> str:
        return self.smectics_dir + "\\smectic_params_" + self.get_density_from_path(source_path) + ".txt"


class CorrelationLengthAnalyzer(ParameterAnalyzer):
    def __init__(self):
        super().__init__()


    def print_correlation_params_for_all_files(self, smectic_file_prefix: str) -> None:  
        files_to_analyze = self.get_all_files_in_dir_to_analyze(self.smectics_dir, smectic_file_prefix)
        for file in files_to_analyze:
            self.print_correlation_params_for_single_file(file)


    def print_correlation_params_for_single_file(self, file: str) -> None:
        fit_params = self.read_and_fit_data(file)
        
        results = self.prepare_results_to_print(fit_params)
        target = self.target_location(file)
        self.print_results_to_file(results, target)


    def read_smectic_params_data(self, file) -> np.array:
        with open(file, 'r') as f:
            data = np.array([line.split() for line in f], dtype=np.float16)
        
        if data.shape[-1] == 2:
            return data
        else:
            raise ValueError(f"Row is expected to have 2 colums. Current number of rows: {data.shape[-1]}")
        
        
    def read_and_fit_data(self, file) -> list[float, float, float, float]:
        data = self.read_smectic_params_data(file)
        fit_params = CorrelationLength().fit_decay(data)
        return fit_params


    def prepare_results_to_print(self, parameters: list[float]) -> list[tuple[str, str]]:
        keys = ["Amplitude:", "DeltaAmplitude:", "CorrelationLength:", "DeltaCorrelationLength:"]
        str_parameters = [str(param) for param in parameters]
        results = zip(keys, str_parameters)
        return results


    def target_location(self, source_path) -> str:
        return self.correlation_dir + "\\correlation_params_" + self.get_density_from_path(source_path) + ".txt"


class EllipsisAnalyzer(ParameterAnalyzer):
    def __init__(self, N_SLICES: int, N_PIX: int = 150):
        super().__init__()
        self.N_SLICES = N_SLICES
        self.N = N_PIX


    def print_semiaxes_for_all_files(self, director_file_prefix: str) -> None:
        files_to_analyze = self.get_all_files_in_dir_to_analyze(self.directors_dir, director_file_prefix)

        for file in files_to_analyze:
            self.print_semiaxes_for_single_file(file)


    def print_semiaxes_for_single_file(self, file: str) -> None:
        matrix = self.read_directors_data(file)
        parameters = self.calculate_semiaxes_from_matrix(matrix)
        results = self.prepare_results_to_print(parameters)
        target = self.target_location(file)
        self.print_results_to_file(results, target)


    def read_directors_data(self, file) -> DirectorsMatrix:
        matrix = DirectorsMatrix(size=(self.N, self.N))
        matrix.read_matrix_from_file(file)
        return matrix


    def calculate_semiaxes_from_matrix(self, matrix: DirectorsMatrix) -> list[float]:
        screen_width = matrix.size[1]
        step_size = screen_width // self.N_SLICES
        semiaxes = []

        for i in range(2, screen_width-2, step_size):
            slice = matrix.get_slice(x_lim=(i, i + step_size))
            fit_result = self.calculate_semiaxes_for_slice(slice)
            semiaxes.append(fit_result)
        return semiaxes


    def calculate_semiaxes_for_slice(self, slice) -> float:
        return EllipsisSemiaxes().calculate_semiaxes(slice)


    def prepare_results_to_print(self, parameters: list[tuple[float, float]]) -> list[tuple[str, str, str]]:
        slice_positions = [str((i+1/2)/self.N_SLICES) for i in range(self.N_SLICES)]
        ex = [str(param[0]) for param in parameters]
        ey = [str(param[1]) for param in parameters]
        results = zip(slice_positions, ex, ey)
        return results


    def target_location(self, source_path) -> str:
        return self.ellipsis_dir + "\\ellipsis_semiaxes_" + self.get_density_from_path(source_path) + ".txt"
