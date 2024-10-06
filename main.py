from common.PROCESSING import *
from common.TIMERS import *

class GlobalParameters:
    user = "Szymek"
    location = "G:/lammps dane/two_domains/all_snapshots_0.32.lammpstrj"
    density = location.split('_')[-1].split('.')[0] + '.' + location.split('_')[-1].split('.')[1],
    mode = location.split('/')[-2],
    params = {
        'LOCATION' : location,
        'DENSITY' : density,
        'MODE' : mode,
        'TARGET_FILE' : "C:/Users/" + user + "/Desktop/LAMMPS_matrices/directors_screen_bulk_" + mode + '_' + density + ".txt",

        'NP' : 10,
        'BATCH_START' : 3,
        'BATCH_STOP' : 10,
        'DIRECTOR_PERIODS' : 2,

        'ATOMS_IN_MOLECULE' : 11,
        'ANALYSE_WALL' : False,
        'ELIMINATE_GOLDSTONE' : False,
        'PLANE' : "xz",
        'SIZE' : [150, 150],
        'PIXEL_TYPE' : DirectorPixel
    }

    def __init__(self, **kwargs) -> None:
        self.params.update(kwargs)


#TODO: move timer to analyzer
# is PIXEL_TYPE necessary?

if __name__ == '__main__':
    params = GlobalParameters()
    analyzer = DirectorFullAnalyzer(params['LOCATION'])

    screen = Screen(params['size'], params['PIXEL_TYPE'])
    timer_analysis = Timer()
    timer_printer = Timer(message_start='Here comes the heatmap!')

    timer_analysis.print_duration(
        analyzer.get_results_in_parallel(screen, analyzer, params)
    )

    timer_printer.print_duration(
        ScreenPrinter(screen).print_screen(params['TARGET_FILE'])
    )
