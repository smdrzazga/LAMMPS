import banana_lib as sz

location = "G:/lammps dane/two_domains/all_snapshots_0.32.lammpstrj"
target = "G:/lammps dane/two_domains/pokolorowane.lammpstrj"

axis = 'y'


i = 0
d = {'x': 0, 'y':1, 'z':2}
with open(location, "r") as f:
    with open(target, "w") as t:
        for line in f:
            i += 1

            try:
                atom = sz.Atom(line.split()[0], *line.split()[-3:])
            except:
                if "id" in line: print(f"{line[:-1]} R G B", file=t)
                else: print(line, end='', file=t)

                continue
            

            # create molecule every 11 atoms
            if atom.id % 11 == 1:
                molecule = sz.Molecule(atom.id//11 + 1, 11)
            
            molecule.add(atom)

            # if molecule is fully read then
            if atom.id % 11 == 0:
                director = molecule.director()
                for element in molecule.comp:
                    R, G, B = 0, 0, 0
                    if director[d[axis]] < 0: R = -1*director[d[axis]]
                    else: G = director[d[axis]]
                    print(f"{element.id} {element.type} {element.position[0]} {element.position[1]} {element.position[2]} {R:.3f} {G:.3f} {B:.3f}", file=t)

            if i > 2e7:
                break