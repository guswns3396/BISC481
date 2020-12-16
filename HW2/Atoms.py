import numpy as np
from numpy.linalg import norm

class Atoms(object):
    def __init__(self):
        self.children = []
        coords = np.array([
            [26.626,   5.330,   3.074],
            [25.816,   5.062,   4.180],
            [25.813,   5.924,   5.279],
            [26.617,   7.065,   5.270],
            [27.429,   7.329,   4.166],
            [27.432,   6.468,   3.066]
        ])
        coords -= np.mean(coords, axis=0)
        
        theta = np.linspace(0, 3*np.pi, 6, endpoint=False)
        dx = 0.5*np.cos(theta)
        dy = 0.5*np.sin(theta)
        dz = 1.5*dx*dy

        coords[:,0] += dx
        coords[:,1] += dy
        coords[:,2] += dz

        self.coords = coords
    
    def __len__(self):
        return len(self.coords)
    
    def save_pdb(self, fname):
        template = "ATOM  {:>5d}  {:4s}{:3s} A{:4d}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00 99.99           C  \n"
        
        FH = open("{}.pdb".format(fname), "w")
        for i in range(len(self)):
            FH.write(template.format(i+1, "C{:d}".format(i+1), "BNZ", 1, self.coords[i,0], self.coords[i,1], self.coords[i,2]))
        
        FH.write("CONECT{:5d}{:5d}{:5d}\n".format(1, 2, 6))
        FH.write("CONECT{:5d}{:5d}{:5d}\n".format(2, 3, 1))
        FH.write("CONECT{:5d}{:5d}{:5d}\n".format(3, 4, 2))
        FH.write("CONECT{:5d}{:5d}{:5d}\n".format(4, 5, 3))
        FH.write("CONECT{:5d}{:5d}{:5d}\n".format(5, 6, 4))
        FH.write("CONECT{:5d}{:5d}{:5d}\n".format(6, 1, 5))
        FH.close()
        

atoms = Atoms()

