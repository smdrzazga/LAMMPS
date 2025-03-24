import banana_lib as sz
import numpy as np
from matplotlib import cm
from matplotlib.colors import hsv_to_rgb


def palette(director, polarization):
    director_proj = director[[0,1]]
    polarization_proj = polarization[[0,1]]
    value = np.cross(director_proj, polarization_proj)

    R, G, B = 0, 0, 0
    if value > 0:
        R = value
    else:
        G = abs(value)

    return R, G, B


# location = r"G:\lammps dane\two_domains\all_snapshots_1.01.lammpstrj"
# target = r"G:\lammps dane\two_domains\coloured_domains_1.01.lammpstrj"
location = r"C:\Users\Szymek\Desktop\middle_snapshot_22500000.lammpstrj"
target = r"C:\Users\Szymek\Desktop\coloured_domains_1.01.lammpstrj"


num_bananas = 528000

i = 0
d = {'x': 0, 'y':1, 'z':2}
with open(location, "r") as f:
    with open(target, "w") as t:
        for line in f:
            i += 1

            try:
                l = line.split()
                atom = sz.Atom(l[0], l[-3:], type=l[1])
            except:
                if "id" in line: print(f"{line[:-1]} R G B", file=t)
                else: print(line, end='', file=t)

                continue
            
            if atom.id <= num_bananas:
                # create molecule every 11 atoms
                if atom.id % 11 == 1:
                    molecule = sz.Molecule(atom.id//11 + 1, 11)
                
                molecule.add(atom)

                # if molecule is fully read then
                if atom.id % 11 == 0:
                    director = molecule.director()
                    polarization = molecule.polarization()
                    for element in molecule.comp:
                        R, G, B = palette(director, polarization)
                        # print(f"{element.id} {molecule.id} {element.type} {element.position[0]} {element.position[1]} {element.position[2]} {R:.3f} {G:.3f} {B:.3f}", file=t)
                        print(f"{element.id} {element.type} {element.position[0]} {element.position[1]} {element.position[2]} {R:.3f} {G:.3f} {B:.3f}", file=t)

            # if i > 1e7:
            #     break
