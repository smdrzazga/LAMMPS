import banana_lib as sz
import numpy as np

def palette(director):
    R = director[0]
    B = -0.5*director[0] - 1.73/2*director[1]
    G = -0.5*director[0] + 1.73/2*director[1]    

    return R, G, B

location = r"E:\LAMMPS_data\two_domains\peter_case\all_snapshots_0.32.lammpstrj"
target = r"E:\LAMMPS_data\two_domains\peter_case\b.lammpstrj"
# location = r"C:\Users\komok\Desktop\new\all_snapshots_0.32.lammpstrj"
# target = r"C:\Users\komok\Desktop\new\a.lammpstrj"
# location = "C:/Users/komok/Desktop/middle_snapshot_9500000.lammpstrj"
# target = "C:/Users/komok/Desktop/a.lammpstrj"

axis = 'x'
num_bananas = 264000

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
                        R, G, B = palette(director)
                        print(f"{element.id} {element.type} {element.position[0]} {element.position[1]} {element.position[2]} {R:.3f} {G:.3f} {B:.3f}", file=t)

            else:
                if atom.id % 10 == 1:
                    molecule = sz.Molecule(atom.id//10 + 1, 10)
                
                molecule.add(atom)

                # if molecule is fully read then
                if atom.id % 10 == 0:
                    director = molecule.director()
                    for element in molecule.comp:
                        R, G, B = palette(director)
                        print(f"{element.id} {element.type} {element.position[0]} {element.position[1]} {element.position[2]} {R:.3f} {G:.3f} {B:.3f}", file=t)

            # if i > 35e6:
            #     break

            # SPLIT INTO 10 AND 11 ATOM MOLECULES