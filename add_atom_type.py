location = r'C:\Users\Szymek\Desktop\two_domains\tworzenie_atomow_z_pliku.txt'
target = r'C:\Users\Szymek\Desktop\two_domains_start.txt'

with open(location, "r+") as f:
    with open(target, "w+") as t:            
        for line in f:
            l = line.split()
            try:
                id = int(l[0])
                type = int(l[2])
            except:
                print(line, end='', file=t)
                continue

            if id > 66000:
                print(f"{id} {l[1]} {type+1} {' '.join(l[3:])}", file=t)
            else:
                print(line, end='', file=t)
