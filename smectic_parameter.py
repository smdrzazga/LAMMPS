import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from multiprocessing.pool import Pool
import mmap
import time
from scipy.optimize import curve_fit


locations = [
    # "E:/LAMMPS_data/6k/all_snapshots_0.318.lammpstrj",    
    # "E:/LAMMPS_data/6k/all_snapshots_0.316.lammpstrj",
    # "E:/LAMMPS_data/6k/all_snapshots_0.315.lammpstrj",
    # "E:/LAMMPS_data/6k/all_snapshots_0.314.lammpstrj",
    # "E:/LAMMPS_data/6k/all_snapshots_0.312.lammpstrj",
    # "E:/LAMMPS_data/6k/all_snapshots_0.31.lammpstrj",
    # "E:/LAMMPS_data/6k/all_snapshots_0.308.lammpstrj",
    # "E:/LAMMPS_data/6k/all_snapshots_0.305.lammpstrj",
    "E:/LAMMPS_data/6k/all_snapshots_0.3.lammpstrj",
    "E:/LAMMPS_data/6k/all_snapshots_0.29.lammpstrj",
    "E:/LAMMPS_data/6k/all_snapshots_0.28.lammpstrj",
    "E:/LAMMPS_data/6k/all_snapshots_0.32.lammpstrj"
]
target = "C:/Users/komok/Desktop/smectic_params_new.txt"
# target = "C:/Users/komok/Desktop/xd.txt"

NP = 6
SIZE = mmap.ALLOCATIONGRANULARITY * 1000

# input bin size and side of simulation box
SLIZE_THICCNESS = 3
x = 150
z = 150
plane = "xz"
DIRECTOR_PERIODS = 1
SMECTIC_PERIODS = 4


def analyze_batch(n, location):
    screen = sz.Screen(x, z, sz.CenterPixel)
    screenshot = sz.Screenshot(x, z, sz.CenterPixel)
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

    molecule_len = 11
    C_left = 0 + 0j
    C_right = 0 + 0j
    C = 0 + 0j

    with open(location, "r+") as f:
        data = mmap.mmap(f.fileno(), length = SIZE, offset = (n)*SIZE)
        molecule = sz.Molecule(1, 11)        
        while True:
            line = data.readline().decode()

            if not line:
                break

            # read atoms one by one from file 
            try:
                atom = sz.Atom( line.split()[0], *line.split()[-3:])
            except:
                continue

            # create molecule every 11 atoms
            if atom.id % 11 == 1:
                molecule = sz.Molecule(atom.id//11 + 1, molecule_len)
            
            molecule.add(atom)

            # if molecule is fully read then
            if len(molecule.comp) == molecule.atoms:
                # translate molecule center back to simulation box 
                center = sz.Atom(atom.id, *molecule.center_of_mass())
                center = sz.wrap_atom_to_box(center, box)
        
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)
                C += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                if center.position[0] < box.x//2:
                    C_left += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                else:
                    C_right += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                
                # reject if center is not close to the wall, else add to screenshot
                # if molecule.center_of_mass()[0] < 5:
                screenshot.assign(center, *screenshot.determine_pixel(center, box, plane))    


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                # flow_left = box.z / DIRECTOR_PERIODS * np.angle(C_left) / (2*np.pi)
                # flow_right = box.z / DIRECTOR_PERIODS * np.angle(C_right) / (2*np.pi)
                flow = box.z / DIRECTOR_PERIODS * np.angle(C) / (2*np.pi)
                # pix_to_scroll_left = screenshot.pixels_to_scroll(z, box, flow_left)
                # pix_to_scroll_right = screenshot.pixels_to_scroll(z, box, flow_right)
                pix_to_scroll = screenshot.pixels_to_scroll(z, box, flow)      
                # screenshot.scroll(pix_to_scroll_left, side="l")
                # screenshot.scroll(pix_to_scroll_right, side="r")
                screenshot.scroll(pix_to_scroll, side="both")
                
                # add corrected screenshot to the final image
                screen.append_screenshot(screenshot)

                # clear current screenshot and number C measuring flow 
                screenshot = sz.Screenshot(x, z, sz.CenterPixel)
                C_left = 0 + 0j
                C_right = 0 + 0j
                C = 0 + 0j

                # print(pix_to_scroll_left, pix_to_scroll_right)
                # print(pix_to_scroll_both)

    print(f"Task is done!: {n}")
    return screen

def create_plot(screen: sz.Screen, location):
    # instantiate box for getting its size
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

    # get packing fraction from file name 
    density = location.split("_")[-1].split(".")
    density = f"{density[0]}.{density[1]}"

    # create lists for plotting local smectic parameters
    smectic_params = [sz.SmecticParameter(SMECTIC_PERIODS) for i in range(x//SLIZE_THICCNESS)]
    # smectic_params = [sz.SmecticParameter(SMECTIC_PERIODS) for i in range(3)]
    
    for n in range(len(smectic_params)):
        smectic_params[n].read_screen(screen, n*SLIZE_THICCNESS, (n+1)*SLIZE_THICCNESS)
        # print(f"{smectic_params[n].parameter} | {smectic_params[n].count}")
        smectic_params[n].normalize()

    # data for plotting 
    # print(smectic_params)
    range_start = 1
    range_stop = len(smectic_params) - 1
    y_axis = [np.abs(smectic_params[i].parameter) for i in range(range_start, range_stop)]
    x_axis = [i*SLIZE_THICCNESS * 0.1*box.x / 150 for i in range(range_start, range_stop//2)]
    x_axis.extend([(range_stop - i)*SLIZE_THICCNESS * 0.1*box.x / 150 for i in range(range_stop//2, range_stop)])

    # curve fitting
    popt, pvar = curve_fit(lambda x, A, L: A*np.exp(-x / L), x_axis, y_axis, sigma=[0.005]*len(x_axis), absolute_sigma=True)
    popt[1] *= box.x / x     # rescale pixels to atom diameters
    print(f"A = {popt[0]}, \nlambda = {popt[1]}")
    with open(target, "a") as t:
        print(f"d={density}", file=t)
        print(f"A = {popt[0]}, \nvarA = {np.sqrt(pvar[0][0])}, \nlambda = {popt[1]} \nvarLambda = {np.sqrt(pvar[1][1])} \n", file=t)

    # exponential plot
    plt.figure()
    plt.scatter(x_axis , y_axis)
    plt.xlabel(f"distance to wall  [molecule lengths]")
    plt.ylabel(f"tau")
    plt.title(f"smectic order parameter, d={density}")
    plt.savefig(f"smectic_order_parameter_{density}.jpg")
    # plt.show()



if __name__ == '__main__':

    for location in locations:
        screen = sz.Screen(x, z, sz.CenterPixel)
        n = 90
        t1 = time.time()
        with Pool(NP) as executor:
            for result in executor.starmap(analyze_batch, zip([i for i in range(n)], [location]*n)):
            # for result in executor.starmap(analyze_batch, [(locations[0], i) for i in range(n)]):
                screen.append_screenshot(result)
        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")
        print("Here comes the heatmap!")
        create_plot(screen, location)
        print("Bye bye, heatmap!")

