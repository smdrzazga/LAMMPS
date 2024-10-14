import numpy as np

location = r'C:\Users\Szymek\Desktop\wca_cut\WCA_cut.txt'

class Potential:
    def calcEnergy(self, r):
        return None
    
    def calcForce(self, r):
        return None
    
    def point(self, i: int, r: float):
        return f"{i} {r:.3f} {self.calcEnergy(r):.5f} {self.calcForce(r):.5f}"

class LinearPotential(Potential):
    def calcEnergy(self, r):
        return 1000*(1 - r)
    
    def calcForce(self, r):
        return 1000

class WCAPotential(Potential):
    def calcEnergy(self, r):
        return 1/r**12 - 1/r**6
    
    def calcForce(self, r):
        return -6*(r**6 - 2) / (r**13)

N = 5000
r = 0
step = 0.001

header = f"""\
# DATE: 2024-10-14  UNITS: lj\n
WCA_cut
N {N} R {step} {N*step}\n
"""

with open(location, 'w+') as f:
    print(header, file=f)
    for i in range(1, N+1):
        r += step
        if r < 1:   
            print(LinearPotential().point(i, r), file=f)
        else:
            print(WCAPotential().point(i, r), file=f)

print("All done!")