from common.PROCESSING import *
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import time


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


#TODO: use GlobalParameters
if __name__ == '__main__':
    params = GlobalParameters()
    screen = Screen(params['size'], DirectorPixel)

    t1 = time.time()
    with Pool(NP) as executor:
        IDs = [ID for ID in range(NP)]
        analyze = BatchAnalyzer(location).analyze_batch
        for result in executor.starmap(analyze, IDs):
            screen.append_screenshot(result)


        print("Here comes the heatmap!")
        ScreenPrinter(screen).print_screen(target_file)
        print("Finished! Yay!")

    t2 = time.time()
    print(f"Time elapsed: {t2 - t1}")
    print("Bye bye, heatmap!")
