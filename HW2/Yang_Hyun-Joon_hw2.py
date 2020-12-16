# These statements import various python modules which you will not be using 
# directly.
import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize

# This is a custom module created for this HW. Please be sure that 'Atoms.py' is
# in the same directory as this script when running. "atoms" is an object that 
# stores the coordinates for some atoms that make up a molecule we will be minimizing
from Atoms import atoms

PI = 3.1415926
# Here are some parameters in our highly simplified force-field. These will be used in the functions
# below to compute energies
L0 = 1.4
TH0 = 2*PI/3
R0 = 1.4
DW = 1/6

################################################################################
# The following functions define various contributions to the potential energy.
# Please try to understand what they are doing, but you don't need to 
# understand every step! Don't edit this code at all!

def V_length(atoms):
    """This function takes the group of atoms and computes the total bond-length
    contribution to the potential energy"""
    
    Vl = 0 # this is the variable we will store the sum of all the energies in
    N = len(atoms)
    for i in range(N):
        j = (i+1) % N
        length = norm(atoms.coords[i] - atoms.coords[j]) # norm computes the length of a vector
        
        Vl += (length - L0)**2
    
    return Vl

def V_angles(atoms):
    """This function takes the group of atoms and computes the total bond-angle
    contribution to the potential energy"""
    
    Va = 0 # this is the variable we will store the sum of all the energies in
    N = len(atoms)
    for i in range(len(atoms)):
        j = (i+1) % N
        k = (i-1) % N
        x_ij = atoms.coords[j] - atoms.coords[i] # vector from atom i to j
        x_ik = atoms.coords[k] - atoms.coords[i] # vector from atom i to k
        theta = np.arccos(np.dot(x_ij, x_ik)/(norm(x_ij)*norm(x_ik))) # angle between the above two
        
        Va += (theta - TH0)**2
    
    return Va

def V_bonds(atoms):
    """This function takes the group of atoms and computes the total bond energy
    based on bond length and bond angles"""
    
    Vb = V_length(atoms) + V_angles(atoms)
    
    return Vb

def V_lennard_jones(atoms):
    """This function takes the group of atoms and computes the total 
    Lennard-Jones energy"""
    
    Vw = 0 # this is the variable we will store the sum of all the energies in
    N = len(atoms)
    for i in range(N):
        for j in range(i+1, N):
            r = norm(atoms.coords[i] - atoms.coords[j]) # distance from atom i to atom j
            
            Vw += DW*((R0/r)**12 -2*(R0/r)**6) # the Lennard-Jones interaction!
    
    return Vw

################################################################################
# The following function minimizes the given energy using a numerical 
# minimization approach. Don't edit this code at all!

def minimize_V(atoms):
    x0 = atoms.coords.flatten()
    
    def V(x):
        atoms.coords = x.reshape(6,3)
        return V_total(atoms)
    
    res = minimize(V, x0, method='Nelder-Mead', tol=1e-6)
    atoms.coords = res.x.reshape(6, 3)


################################################################################
# YOUR WORK BEGINS HERE!

# STEP 1 - complete this function to compute the total potential energy. It takes the object "atoms"
# as input and returns a number.
def V_total(atoms):
    
    Vt = V_bonds(atoms) + V_lennard_jones(atoms)
    
    return Vt

# STEP 2 - using the "print" function, compute the total energy for the group
# of atoms in their initial state, and report the value you got.

print('Initial energy:', V_total(atoms))

# Now let's save the initial configuration of the atoms using a method of the 
# atoms object. This will create a file called "initial.pdb"
atoms.save_pdb("initial")

# Now we are going to minimize the potential energy, V_total
minimize_V(atoms)

# STEP 3 - print the new minimized energy, and save the atom configuration like 
# before, calling it "final.pdb"
print('Final energy:', V_total(atoms))
atoms.save_pdb("final")

# Hyun-Joon Yang
