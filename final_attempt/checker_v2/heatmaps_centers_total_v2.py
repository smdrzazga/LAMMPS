import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import mmap
import time


location = "E:/LAMMPS_data/6k/all_snapshots_0.32.lammpstrj"


size = mmap.ALLOCATIONGRANULARITY * 200
# input bin size and side of simulation box
x = 150
z = 150
plane = "xz"
periods = 1

colorbar_limits = (500, 1000)
box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

def analyze_batch(n):
    screen = sz.Screen(x, z, sz.CenterPixel)
    screenshot = sz.Screenshot(x, z, sz.CenterPixel)

    molecule = sz.Molecule(1, 11)
    C_left = 0 + 0j
    C_right = 0 + 0j
    C = 0 + 0j
    with open(location, "r+") as f:
        data = mmap.mmap(f.fileno(), length = size, offset = (n)*size)
        molecule = sz.Molecule(1, 11)
        while True:
            # atom = sz.Atom( *read_atom(line).T)
            line = data.readline().decode()

            if not line:
                break

            # read atoms one by one from file 
            try:
                # atom = sz.Atom( *line.strip().split() )
                atom = sz.Atom( line.split()[0], *line.split()[-3:])
            except:
                continue


            # create molecule every 11 atoms
            if atom.id % 11 == 1:
                molecule = sz.Molecule(atom.id//11 + 1, 11)
            
            molecule.add(atom)

            # if molecule is fully read then
            if len(molecule.comp) == molecule.atoms:
                # translate molecule center back to simulation box 
                # center = sz.Atom(atom.id, 'A', *molecule.center_of_mass())
                center = sz.Atom(atom.id, *molecule.center_of_mass())
                center = sz.wrap_atom_to_box(center, box)
        
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)\
                C += molecule.polarization()[1] * np.exp(2j*periods*np.pi * center.position[2] / box.z)
                if center.position[0] < box.x//2:
                    C_left += molecule.polarization()[1] * np.exp(2j*periods*np.pi * center.position[2] / box.z)
                else:
                    C_right += molecule.polarization()[1] * np.exp(2j*periods*np.pi * center.position[2] / box.z)

                # reject if center is not close to the wall, else add to screenshot
                # if molecule.center_of_mass()[0] < 5:
                screenshot.assign(center, *screenshot.determine_pixel(center, box, plane))    


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                flow_left = box.z / periods * np.angle(C_left) / (2*np.pi)
                flow_right = box.z / periods * np.angle(C_right) / (2*np.pi)
                pix_to_scroll_left = screenshot.pixels_to_scroll(z, box, flow_left)
                pix_to_scroll_right = screenshot.pixels_to_scroll(z, box, flow_right)
                # screenshot.scroll(pix_to_scroll_left, side="l")
                # screenshot.scroll(pix_to_scroll_right, side="r")

                flow = box.z / periods * np.angle(C) / (2*np.pi)
                pix_to_scroll_both = screenshot.pixels_to_scroll(z, box, flow)
                screenshot.scroll(pix_to_scroll_both, side="both")

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshot)

                # clear current screenshot and number C measuring flow 
                screenshot = sz.Screenshot(x, z, sz.CenterPixel)
                C_left = 0 + 0j
                C_right = 0 + 0j
                C = 0 + 0j

                print(pix_to_scroll_left, pix_to_scroll_right)
                # print(pix_to_scroll_both)

    print(f"Task is done!: {n}")
    return screen


def create_heatmap(screen: sz.Screen, num_batches: int):
    plt.figure()
    plt.imshow(screen.colour()[:, 4:-4], interpolation="bilinear", clim = (None, None), extent=(0, box.x, 0, box.z), aspect="auto")
    # if iter == 0:
    plt.colorbar()

    plt.xlabel(f"{plane[0]}  [atom diameters]")
    plt.ylabel(f"{plane[1]}  [atom diameters]")

    density = location.split("_")[-1].split(".")
    density = density[0] + "." + density[1]
    step = location.split('_')[-1].split('.')[0]
    plt.title(f"centers, density={density}")

    # show OR save plot - not both at once!
    plt.show()
    # plt.savefig(f"centers_{density}_xz_total.png")
    # plt.savefig(f"centers_starting_conf_1x.png")
    # plt.savefig(f"centers_step_{step}_1x.png")
    




if __name__ == '__main__':
    # locations = [
    #     "C:/Users/Szymek/Desktop/all_snapshots_0.32.lammpstrj"
    # ]

    screen = sz.Screen(x, z, sz.CenterPixel)
    n = 50
    t1 = time.time()
    with ProcessPoolExecutor(6) as executor:
        for result in executor.map(analyze_batch, [i for i in range(n)]):
        # for result in executor.starmap(analyze_batch, [(locations[0], i) for i in range(n)]):
            screen.append_screenshot(result)
        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")
        print("Here comes the heatmap!")
        create_heatmap(screen, n)
        print("Bye bye, heatmap!")

