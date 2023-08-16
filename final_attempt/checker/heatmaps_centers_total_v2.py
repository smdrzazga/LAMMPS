import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import mmap
import time

# locations = [
#     "/lu/topola/home/iq2djt3m/ntb_density_decrease/density_decrease/scaling_decrease/bulk_ntb_start/6k_molecules/step_0.001/d_0.318/all_snapshots_0.318.lammpstrj",
#     "/lu/topola/home/iq2djt3m/ntb_density_decrease/density_decrease/scaling_decrease/bulk_ntb_start/6k_molecules/step_0.001/d_0.317/all_snapshots_0.317.lammpstrj",
#     "/lu/topola/home/iq2djt3m/ntb_density_decrease/density_decrease/scaling_decrease/bulk_ntb_start/6k_molecules/step_0.001/d_0.316/all_snapshots_0.316.lammpstrj",
#     "/lu/topola/home/iq2djt3m/ntb_density_decrease/density_decrease/scaling_decrease/bulk_ntb_start/6k_molecules/step_0.001/d_0.309/all_snapshots_0.309.lammpstrj",
#     "/lu/topola/home/iq2djt3m/ntb_density_decrease/density_decrease/scaling_decrease/bulk_ntb_start/6k_molecules/step_0.001/d_0.308/all_snapshots_0.308.lammpstrj",
#     "/lu/topola/home/iq2djt3m/ntb_density_decrease/density_decrease/scaling_decrease/bulk_ntb_start/6k_molecules/step_0.001/d_0.307/all_snapshots_0.307.lammpstrj"
#     ]

size = mmap.ALLOCATIONGRANULARITY * 200
# input bin size and side of simulation box
x = 150
z = 150
plane = "xz"
periods = 4

colorbar_limits = (1200, 1600)

# location = "C:/Users/Szymek/Desktop/all_snapshots_0.32.lammpstrj"
location = "G:/lammps dane/12k_all_snapshots_0.308.lammpstrj"
def analyze_batch(n):
    screen = sz.Screen(x, z, sz.CenterPixel)
    screenshot = sz.Screenshot(x, z, sz.CenterPixel)
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

    i = 0
    C = 0 + 0j
    with open(location, "r+") as f:
        data = mmap.mmap(f.fileno(), length = size, offset = (n+700)*size)
        molecule = sz.Molecule(1, 11)
        while True:
            line = data.readline().decode()
            # print(line)

            if not line:
                break

            # read atoms one by one from file 
            try:
                atom = sz.Atom( *line.strip().split() )
                i += 1
            except:
                continue

            # create molecule every 11 atoms
            if atom.id % 11 == 1:
                molecule = sz.Molecule(atom.id//11 + 1, 11)
            
            molecule.add(atom)

            # if molecule is fully read then
            if len(molecule.comp) == molecule.atoms:
                # translate molecule center back to simulation box 
                center = sz.Atom(atom.id, 'A', *molecule.center_of_mass())
                center = sz.wrap_atom_to_box(center, box)
        
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)
                C += molecule.polarization()[1] * np.exp(2j*periods*np.pi * center.position[2] / box.z)
                
                # reject if center is not close to the wall, else add to screenshot
                # if molecule.center_of_mass()[0] < 5:
                screenshot.assign(center, *screenshot.determine_pixel(center, box, plane))    


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                flow = box.z / periods * np.angle(C) / (2*np.pi)
                pix_to_scroll = screenshot.pixels_to_scroll(z, box, flow)
                screenshot.scroll(pix_to_scroll)

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshot)

                # clear current screenshot and number C measuring flow 
                screenshot = sz.Screenshot(x, z, sz.CenterPixel)
                C = 0

                print(pix_to_scroll)

    print(f"Task is done!: {n}")
    return screen

def create_heatmap(screen):
    # plt.imshow(screen.colour(), interpolation="none", clim = colorbar_limits)
    plt.figure()
    plt.imshow(screen.colour()[:, 4:-4])
    # if iter == 0:
    plt.colorbar()

    plt.xlabel(f"{plane[0]}")
    plt.ylabel(f"{plane[1]}")

    density = location.split("_")[-1].split(".")
    density = density[0] + "." + density[1]
    plt.title(f"centers, total, d={density}")

    # show OR save plot - not both at once!
    # plt.show()
    plt.savefig(f"centers_{density}_xz_total.png")




if __name__ == '__main__':
    # locations = [
    #     "C:/Users/Szymek/Desktop/all_snapshots_0.32.lammpstrj"
    # ]

    screen = sz.Screen(x, z, sz.CenterPixel)
    n = 200
    t1 = time.time()
    with ProcessPoolExecutor(10) as executor:
        for result in executor.map(analyze_batch, [i for i in range(n)]):
        # for result in executor.starmap(analyze_batch, [(locations[0], i) for i in range(n)]):
            screen.append_screenshot(result)
        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")
        print("Here comes the heatmap!")
        create_heatmap(screen)
        print("Bye bye, heatmap!")

