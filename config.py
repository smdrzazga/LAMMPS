class ProcessingParameters:
        params = {
        'NP' : 10,
        'ANALYZE_RANGE' : (0., 1.)
        }

class FileParameters:
    _user = "Szymek"
    _location = "G:/lammps dane/two_domains/all_snapshots_0.32.lammpstrj"
    _density = _location.split('_')[-1].split('.')[0] + '.' + _location.split('_')[-1].split('.')[1]
    _mode = _location.split('/')[-2]

    params = {
        'LOCATION' : _location,
        'DENSITY' : _density,
        'MODE' : _mode,
        'TARGET_FILE' : "C:/Users/" + _user + "/Desktop/LAMMPS_matrices/directors_screen_bulk_" + _mode + '_' + _density + ".txt"
    }   

class ShapshotParams:
    params = {
        'ANALYSE_WALL' : False,
        'ELIMINATE_GOLDSTONE' : False,
        'DIRECTOR_PERIODS' : 2,
        'ATOMS_IN_MOLECULE' : 11,
        'PLANE' : "xz",
        'SIZE' : (150, 150)
    }

class MiscParams:
    params = {}


class GlobalParameters:
    proc_params = ProcessingParameters.params
    file_params = FileParameters.params
    snap_params = ShapshotParams.params
    misc_params = MiscParams.params

    def __init__(self, **kwargs) -> None:
        for key, value in kwargs:
            if key in self.proc_params.keys():
                self.proc_params.update({key : value})
            elif key in self.file_params.keys():
                self.file_params.update({key : value})
            elif key in self.snap_params.keys():
                self.snap_params.update({key : value})
            else:
                self.misc_params.update({key : value})
