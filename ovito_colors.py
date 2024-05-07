import banana_lib as sz
import numpy as np
from matplotlib import cm
from matplotlib.colors import hsv_to_rgb

# def palette(director):
#     R = director[0]
#     B = -0.5*director[0] - 1.73/2*director[1]
#     G = -0.5*director[0] + 1.73/2*director[1]    
#     return R, G, B

def palette(director, colors, N):
    cast = director[0] + 1.j*director[1]
    angle = (np.angle(cast) / 2 / np.pi) + 0.5       # angle ranging from 0 to 1 in full angle units

    i = int(N * angle) - 1
    [R, G, B] = colors[i, :]

    return R, G, B


location = r"G:\lammps dane\4z_2x\all_snapshots_0.305.lammpstrj"
target = r"G:\lammps dane\4z_2x\coloured_4z_2x_0.305.lammpstrj"

axis = 'x'
num_bananas = 528000

N = 255
hsv = np.dstack((np.linspace(0, 1, N), np.ones(N), np.ones(N)))
colors = hsv_to_rgb(hsv.reshape((N, 3)))

i = 0
d = {'x': 0, 'y':1, 'z':2}
with open(location, "r") as f:
    with open(target, "w") as t:
        for line in f:
            i += 1

            try:
                l = line.split()
                atom = sz.Atom(l[0], *l[-3:], type=l[1])
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
                    for element in molecule.comp:
                        R, G, B = palette(director, colors, N)
                        print(f"{element.id} {element.type} {element.position[0]} {element.position[1]} {element.position[2]} {R:.3f} {G:.3f} {B:.3f}", file=t)

            else:
                if atom.id % 10 == 1:
                    molecule = sz.Molecule(atom.id//10 + 1, 10)
                
                molecule.add(atom)

                # if molecule is fully read then
                if atom.id % 10 == 0:
                    director = molecule.director()
                    for element in molecule.comp:
                        R, G, B = palette(director, colors, N)
                        print(f"{element.id} {element.type} {element.position[0]} {element.position[1]} {element.position[2]} {R:.3f} {G:.3f} {B:.3f}", file=t)

            # if i > 1e7:
            #     break

            # SPLIT INTO 10 AND 11 ATOM MOLECULES