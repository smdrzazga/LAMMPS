_location = "G:/lammps dane/6k/all_snapshots_0.32.lammpstrj"

class Parameters:
    params = {}

    def __init__(self, **kwargs) -> None:
        for key, value in kwargs.items():
            if key in self.params.keys():
                self.params.update({key : value})

class ProcessingParameters(Parameters):
        params = {
            'INPUT_FILE' : None,
            'NP' : 10,
            'ANALYZE_RANGE' : (0., 0.005)
        }

class FileParameters(Parameters):
    _user = "Szymek"
    _location = 'C:/Users/Szymek/Desktop/middle_snapshot_4000000.lammpstrj'
    _density = _location.split('_')[-1].split('.')[0] + '.' + _location.split('_')[-1].split('.')[1]
    _mode = _location.split('/')[-2]

    params = {
        'INPUT_FILE' : _location,
        'DENSITY' : _density,
        'MODE' : _mode,
        'OUTPUT_FILE' : "C:/Users/" + _user + "/Desktop/TEST_DIRECTORS" + _mode + '_' + _density + ".txt"
    }   

class ShapshotParams(Parameters):
    params = {
        'ANALYSE_WALL' : False,
        'ELIMINATE_GOLDSTONE' : False,
        'DIRECTOR_PERIODS' : 2,
        'ATOMS_IN_MOLECULE' : 11,
        'PLANE' : "xz",
        'SIZE' : (150, 150)
    }

class MiscParams(Parameters):
    pass


class GlobalParameters:
    proc_params = ProcessingParameters.params
    file_params = FileParameters.params
    snap_params = ShapshotParams.params
    misc_params = MiscParams.params

    def __init__(self, **kwargs) -> None:
        for key, value in kwargs.items():
            if key in self.proc_params.keys():
                self.proc_params.update({key : value})
            if key in self.file_params.keys():
                self.file_params.update({key : value})
            if key in self.snap_params.keys():
                self.snap_params.update({key : value})
            else:
                self.misc_params.update({key : value})
