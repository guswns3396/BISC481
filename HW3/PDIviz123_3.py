"""
The MIT License (MIT)

Copyright (c) 2015 Judemir Ribeiro, Francisco Melo and Andreas Schueller

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from pymol import *
import pymol as pm
import os,glob,re,colorsys,sys
#from string import maketrans
import tkinter
import tkinter as tk
import tkinter.colorchooser
import tkinter.filedialog
import Pmw
from pymol import cmd, cgo, CmdException
import numpy as np

##Fixed atom vdw radii dict
#reference: Chothia 1975 and NACCESS
stored.PDI_Viz_res_dict = {
           'ILE': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'O': 1.4, 'N': 1.65, 'CD1': 1.87, 'CG1': 1.87, 'CG2': 1.87},
           'GLN': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.87, 'O': 1.4, 'N': 1.65, 'CD': 1.76, 'NE2': 1.65, 'OE1': 1.4},
           'GLX': {'AE1': 1.5, 'C': 1.76, 'AE2': 1.5, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'O': 1.4, 'N': 1.65, 'CD': 1.87},
           'GLY': {'OXT': 1.4, 'CA': 1.87, 'C': 1.76, 'O': 1.4, 'N': 1.65},
           '2MG': {'C1*': 1.8, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O6': 1.4, 'C5*': 1.8, 'C2A': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C8': 1.8, 'O2*': 1.4, "O2'": 1.4, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'OXT': 1.4, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2': 1.8, 'C2*': 1.8, 'N1': 1.6, 'N2': 1.6, 'N3': 1.6, 'C4*': 1.8, 'N7': 1.6, 'O5*': 1.4, 'N9': 1.6, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4},
           'HEM': {'CMD': 1.9, 'CMC': 1.9, 'CMB': 1.9, 'CMA': 1.9, 'O1D': 1.35, 'FE': 1.47, 'O1A': 1.35, 'N D': 1.55, 'CAC': 1.9, 'CAB': 1.9, 'CAA': 1.9, 'C1D': 1.78, 'N A': 1.55, 'N C': 1.55, 'N B': 1.55, 'C4D': 1.78, 'CGA': 1.9, 'CGD': 1.9, 'C2B': 1.78, 'C2C': 1.78, 'C2A': 1.78, 'C4A': 1.78, 'C4B': 1.78, 'C4C': 1.78, 'CAD': 1.9, 'CHA': 2.0, 'CHB': 2.0, 'CHC': 2.0, 'CHD': 2.0, 'C2D': 1.78, 'OXT': 1.4, 'O2D': 1.35, 'O2A': 1.35, 'CBB': 1.9, 'CBC': 1.9, 'CBA': 1.9, 'CBD': 1.9, 'C1C': 1.78, 'C1B': 1.78, 'C1A': 1.78, 'C3A': 1.78, 'C3C': 1.78, 'C3B': 1.78, 'C3D': 1.78},
           'GLU': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.87, 'O': 1.4, 'N': 1.65, 'OE2': 1.4, 'CD': 1.76, 'OE1': 1.4},
           'CYS': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'O': 1.4, 'N': 1.65, 'SG': 1.85},
           'ASP': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'O': 1.4, 'N': 1.65, 'OD1': 1.4, 'OD2': 1.4},
           'SER': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'O': 1.4, 'N': 1.65, 'OG': 1.4},
           'LYS': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.87, 'O': 1.4, 'N': 1.65, 'NZ': 1.5, 'CE': 1.87, 'CD': 1.87},
           'PRO': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.87, 'O': 1.4, 'N': 1.65, 'CD': 1.87},
           'ASX': {'C': 1.76, 'AD2': 1.5, 'AD1': 1.5, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'O': 1.4, 'N': 1.65},
           'ASN': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'O': 1.4, 'N': 1.65, 'OD1': 1.4, 'ND2': 1.65},
           'VAL': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'O': 1.4, 'N': 1.65, 'CG1': 1.87, 'CG2': 1.87},
           'THR': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'OG1': 1.4, 'O': 1.4, 'N': 1.65, 'CG2': 1.87},
           'HOH': {'O': 1.4},
           'HIS': {'C': 1.76, 'CE1': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'O': 1.4, 'N': 1.65, 'CD2': 1.76, 'ND1': 1.65, 'NE2': 1.65},
           'TRP': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'CH2': 1.76, 'O': 1.4, 'N': 1.65, 'CZ2': 1.76, 'CE2': 1.76, 'CE3': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CZ3': 1.76, 'NE1': 1.65},
           'H2U': {'C1*': 1.8, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O4': 1.4, 'O2': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'O2*': 1.4, "O2'": 1.4, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'OXT': 1.4, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2': 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'C4*': 1.8, 'O5*': 1.4, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4},
           'PHE': {'C': 1.76, 'CE1': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'O': 1.4, 'N': 1.65, 'CZ': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE2': 1.76},
           'ALA': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'O': 1.4, 'N': 1.65},
           'MET': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.87, 'O': 1.4, 'N': 1.65, 'CE': 1.87, 'SD': 1.85},
           'ACE': {'OXT': 1.4, 'C': 1.76, 'CH3': 1.87, 'O': 1.4},
           'LEU': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.87, 'O': 1.4, 'N': 1.65, 'CD1': 1.87, 'CD2': 1.87},
           'ARG': {'C': 1.76, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.87, 'NE': 1.65, 'O': 1.4, 'N': 1.65, 'CZ': 1.76, 'NH1': 1.65, 'NH2': 1.65, 'CD': 1.87},
           'TYR': {'C': 1.76, 'CE1': 1.76, 'OH': 1.4, 'OXT': 1.4, 'CB': 1.87, 'CA': 1.87, 'CG': 1.76, 'O': 1.4, 'N': 1.65, 'CZ': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE2': 1.76},
           'DG': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O6': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C8': 1.8, 'C2': 1.8, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N2': 1.6, 'N3': 1.6, 'C4*': 1.8, 'N7': 1.6, 'O5*': 1.4, 'N9': 1.6, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4, "OP3": 1.4}, #OP3 added to match NACCESS
           'DC': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O2': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C2': 1.8, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'N4': 1.6, 'C4*': 1.8, 'O5*': 1.4, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4,"OP3": 1.4}, #OP3 added to match NACCESS
           'DA': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C8': 1.8, 'N6': 1.6, 'C2': 1.8, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'C4*': 1.8, 'N7': 1.6, 'O5*': 1.4, 'N9': 1.6, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4,"OP3": 1.4}, #OP3 added to match NACCESS
           'DT': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O4': 1.4, 'O2': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C2': 1.8, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'C4*': 1.8, 'O5*': 1.4, "O5'": 1.4, 'O3*': 1.4, 'C5M': 1.8, "O3'": 1.4, 'C7': 1.8 , "OP3": 1.4},
           'G': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O6': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C8': 1.8, 'O2*': 1.4, "O2'": 1.4,'C2': 1.8, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N2': 1.6, 'N3': 1.6, 'C4*': 1.8, 'N7': 1.6, 'O5*': 1.4, 'N9': 1.6, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4, "OP3": 1.4}, #OP3 added to match NACCESS
           'C': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O2': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C2': 1.8,'O2*': 1.4, "O2'": 1.4, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'N4': 1.6, 'C4*': 1.8, 'O5*': 1.4, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4,"OP3": 1.4}, #OP3 added to match NACCESS
           'A': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'C8': 1.8, 'N6': 1.6,'O2*': 1.4, "O2'": 1.4, 'C2': 1.8, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'C4*': 1.8, 'N7': 1.6, 'O5*': 1.4, 'N9': 1.6, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4,"OP3": 1.4}, #OP3 added to match NACCESS
           'T': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O4': 1.4, 'O2': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4,'O2*': 1.4, "O2'": 1.4 ,'C2': 1.8, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'C4*': 1.8, 'O5*': 1.4, "O5'": 1.4, 'O3*': 1.4, 'C5M': 1.8, "O3'": 1.4, 'C7': 1.8 , "OP3": 1.4},
           'U': {'C1*': 1.8, 'OP1': 1.4, 'OP2': 1.4, 'O1P': 1.4, "C3'": 1.8, "C1'": 1.8, "C5'": 1.8, 'O4': 1.4, 'O2': 1.4, 'C5*': 1.8, 'O4*': 1.4, 'C3*': 1.8, "O4'": 1.4, 'O2*': 1.4, "O2'": 1.4, 'C6': 1.8, 'C5': 1.8, 'C4': 1.8, 'O2P': 1.4, 'P': 1.9, "C4'": 1.8, "C2'": 1.8, 'C2': 1.8, 'C2*': 1.8, 'N1': 1.6, 'N3': 1.6, 'C4*': 1.8, 'O5*': 1.4, "O5'": 1.4, 'O3*': 1.4, "O3'": 1.4, "OP3": 1.4}} #C7,OP3 added to match NACCESS

#DNA backbone polar atoms
stored.PDI_Viz_dna_polar_bb = {'DA':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'DC':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'DG':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'DT':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'U': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'A': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'C': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'G': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'T': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\'']}
#DNA base polar atoms
stored.PDI_Viz_dna_polar_ba =  {'DA':['N1','N3','N6','N7'],
                          'DC':['N3','O2','N4','N1'],
                          'DG':['O6','N1','N2','N3','N7'],
                          'DT':['O2','O4','N3','N1'],
                          'U': ['O2','O4','N3','N1'],
                          'A':['N1','N7','N3','N6'],
                          'C':['N3','O2','N4','N1'],
                          'G':['O6','N3','N7','N2','N1'],
                          'T':['O2','O4','N3','N1']}
#DNA backbone apolar
stored.PDI_Viz_dna_apolar_bb = {'DA':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'DC':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'DG':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'DT':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'U': ['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'A':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'C':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'G':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'T':['C1\'','C2\'','C3\'','C4\'','C5\'']}
#DNA base apolar
stored.PDI_Viz_dna_apolar_ba =   {'DA':['C2','C4','C5','C6','C8','N9'],
                                  'DC':['C2','C4','C5','C6'],
                                  'DG':['C2','C4','C5','C6','C8','N9'],
                                  'DT':['C2','C4','C5','C6','C7','C5M'],
                                  'U': ['C2','C4','C5','C6',],
                                  'A':['C2','C4','C5','C6','C8','N9'],
                                  'C':['C2','C4','C5','C6'],
                                  'G':['C2','C4','C5','C6','C8','N9'],
                                  'T':['C2','C4','C5','C6','C7','C5M'],}

#DNA base H bond acceptors
stored.PDI_Viz_dna_Hbond_A = {'DA':['N1','N7','N3'],
                              'DC':['N3','O2'],
                              'DG':['O6','N3','N7'],
                              'DT':['O2','O4'],
                              'U': ['O2','O4'],
                              'A':['N1','N7','N3'],
                              'C':['N3','O2'],
                              'G':['O6','N3','N7'],
                              'T':['O2','O4']}
#DNA base H bond donors
stored.PDI_Viz_dna_Hbond_D = {'DA':['N6'],
                          'DC':['N4'],
                          'DG':['N2','N1'],
                          'DT':['N3'],
                          'U' :['N3'],
                          'A':['N6'],
                          'C':['N4'],
                          'G':['N2','N1'],
                          'T':['N3']}
#DNA thymine methyl
stored.PDI_Viz_dna_Tmet = {'DA':[],
                       'DC':[],
                       'DG':[],
                       'DT':['C7','C5M'],
                        'U':[],
                        'A':[],
                       'C':[],
                       'G':[],
                       'T':['C7','C5M']}
#DNA base atoms
stored.PDI_Viz_dna_Hn = {'DA':['C2','C8'],
                     'DC':['C5','C6'],
                     'DG':['C8'],
                     'DT':['C6'],
                     'U': ['C5','C6'],
                     'A':['C2','C8'],
                     'C':['C5','C6'],
                     'G':['C8'],
                     'T':['C6']}
#Protein side chain acceptors
stored.PDI_Viz_prot_Hbond_A = {'GLY':[],
                 'ALA':[],
                 'VAL':[],
                 'LEU':[],
                 'ILE':[],
                 'MET':['SD'],
                 'PRO':[],
                 'PHE':[],
                 'TRP':[],
                 'SER':['OG'],
                 'THR':['OG1'],
                 'GLN':['OE1'],
                 'LYS':[],
                 'TYR':['OH'],
                 'ASN':['OD1'],
                 'CYS':[],
                 'GLU':['OE1','OE2'],
                 'ASP':['OD1','OD2'],
                 'ARG':[],
                 'HIS':['ND1'],
                 'HID':['NE2'],
                 'HIE':['ND1'],
                 'HIP':[]}
#Protein side chain donors
stored.PDI_Viz_prot_Hbond_D = {'GLY':[],
                 'ALA':[],
                 'VAL':[],
                 'LEU':[],
                 'ILE':[],
                 'MET':[],
                 'PRO':[],
                 'PHE':[],
                 'TRP':['NE1'],
                 'SER':['OG'],
                 'THR':['OG1'],
                 'GLN':['NE2'],
                 'LYS':['NZ'],
                 'TYR':['OH'],
                 'ASN':['ND2'],
                 'CYS':['SG'],
                 'GLU':[],
                 'ASP':[],
                 'ARG':['NH1','NH2','NE'],
                 'HIS':['NE2'],
                 'HID':['ND1'],
                 'HIE':['NE2'],
                 'HIP':['NE2','ND1']}

stored.PDI_Viz_prot_polar = {'GLY':['N','O','OXT'],
                             'ALA':['N','O','OXT'],
                             'VAL':['N','O','OXT'],
                             'LEU':['N','O','OXT'],
                             'ILE':['N','O','OXT'],
                             'MET':['SD','N','O','OXT'],
                             'PRO':['N','O','OXT'],
                             'PHE':['N','O','OXT'],
                             'TRP':['NE1','N','O','OXT'],
                             'SER':['OG','N','O','OXT'],
                             'THR':['OG1','N','O','OXT'],
                             'GLN':['OE1','NE2','N','O','OXT'],
                             'LYS':['NZ','N','O','OXT'],
                             'TYR':['OH','N','O','OXT'],
                             'ASN':['OD1','ND2','N','O','OXT'],
                             'CYS':['SG','N','O','OXT'],
                             'GLU':['OE1','OE2','N','O','OXT'],
                             'ASP':['OD1','OD2','N','O','OXT'],
                             'ARG':['NH1','NH2','NE','N','O','OXT'],
                             'HIS':['ND1','NE2','N','O','OXT'],
                             'HID':['ND1','NE2','N','O','OXT'],
                             'HIE':['ND1','NE2','N','O','OXT'],
                             'HIP':['ND1','NE2','N','O','OXT']}

stored.PDI_Viz_prot_apolar = {'GLY':['CA','C',],
                 'ALA':['CA','C','CB'],
                 'VAL':['CA','C','CB','CG1','CG2'],
                 'LEU':['CA','C','CB','CG','CD1','CD2'],
                 'ILE':['CA','C','CB','CG1','CG2','CD1'],
                 'MET':['CA','C','CB','CG','CE'],
                 'PRO':['CA','C','CB','CG','CD'],
                 'PHE':['CA','C','CB','CG','CD1','CD2','CE1','CE2','CZ'],
                 'TRP':['CA','C','CB','CG','CD1','CD2','CE2','CZ2','CH2','CZ3','CE3'],
                 'SER':['CA','C','CB'],
                 'THR':['CA','C','CB','CG2'],
                 'GLN':['CA','C','CB','CG','CD'],
                 'LYS':['CA','C','CB','CG','CD','CE'],
                 'TYR':['CA','C','CB','CG','CD1','CD2','CE1','CE2','CZ'],
                 'ASN':['CA','C','CB','CG'],
                 'CYS':['CA','C','CB'],
                 'GLU':['CA','C','CB','CG','CD'],
                 'ASP':['CA','C','CB','CG'],
                 'ARG':['CA','C','CB','CG','CD','CZ'],
                 'HIS':['CA','C','CB','CG','CD2','CE1'],
                 'HID':['CA','C','CB','CG','CD2','CE1'],
                 'HIE':['CA','C','CB','CG','CD2','CE1'],
                 'HIP':['CA','C','CB','CG','CD2','CE1']}

stored.PDI_Viz_dna_majorgroove = {'DA':['C5', 'C6', 'C8', 'N6', 'N7'],
                                  'DC':['C4', 'C5', 'C6', 'N4'],
                                  'DG':['C5', 'C6', 'C8', 'N7', 'O6'],
                                  'DT':['C4', 'C5', 'C6', 'C7', 'O4','C5M'],
                                  'U':['C4', 'C5', 'C6', 'O4'],
                                  'A':['C5', 'C6', 'C8', 'N6', 'N7'],
                                  'C':['C4', 'C5', 'C6', 'N4'],
                                  'G':['C5', 'C6', 'C8', 'N7', 'O6'],
                                  'T':['C4', 'C5', 'C6', 'C7', 'O4','C5M']} # RNA
stored.PDI_Viz_dna_minorgroove = {'DA':['C2', 'C4', 'N1', 'N3'],
                                  'DC':['C2', 'O2'],
                                  'DG':['C2', 'C4', 'N2', 'N3'],
                                  'DT':['C2', 'N3', 'O2'],
                                   'U':['C2', 'N3', 'O2'],
                                  'A':['C2', 'C4', 'N1', 'N3'],
                                  'C':['C2', 'O2'],
                                  'G':['C2', 'C4', 'N2', 'N3'],
                                  'T':['C2', 'N3', 'O2']} # RNA
stored.PDI_Viz_dna_ambi = {'DA':['N9'],
                           'DC':['N1', 'N3'],
                           'DG':['N1', 'N9'],
                           'DT':['N1'],
                            'U':['N1'],
                           'A':['N9'],
                           'C':['N1', 'N3'],
                           'G':['N1', 'N9'],
                           'T':['N1']} # RNA

#protein backbone and bases
stored.PDI_Viz_prot_backbone = ['N', 'CA', 'C', 'O', 'OXT']
stored.PDI_Viz_prot_base = ['CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CZ', 'CZ2',
                                          'CZ3', 'CH2', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ', 'NH1', 'NH2','OD1', 'OD2', 'OG',
                                          'OG1', 'OG2', 'OE1', 'OE2', 'OH', 'SD', 'SG']
#dna backbone and bases
stored.PDI_Viz_dna_backbone = ['C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'','O2\'', 'O3\'', 'O4\'', 'O5\'', 'OP1', 'OP2', 'OP3', 'P']
stored.PDI_Viz_dna_base = ['C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6','C5M']

stored.PDI_Viz_IsGUIOpen = False
stored.PDI_Viz_IsColorSelOpen = False
stored.PDI_Viz_IsSaveOpen = False
stored.PDI_Viz_IsHelpOpen = False

stored.PDI_Viz_dSASAcutoff = 0.0

#scratch variable storage
stored.PDI_Viz_p_bb = 'Protein_bb'
stored.PDI_Viz_p_sc = 'Protein_sc'
stored.PDI_Viz_p_polar = 'Protein_polar'
stored.PDI_Viz_p_apolar = 'Protein_apolar'
stored.PDI_Viz_DNA_bb = 'DNA_bb'
stored.PDI_Viz_DNA_base = 'DNA_base'
stored.PDI_Viz_Prot = 'protein'
stored.PDI_Viz_DNA = 'dna'
stored.PDI_Viz_DNA_majorg = 'dna_ma'
stored.PDI_Viz_DNA_minorg = 'dna_mi'
stored.PDI_Viz_DNA_amb = 'dna_amb'
stored.PDI_Viz_DNA_polar = 'dna_polar'
stored.PDI_Viz_DNA_apolar = 'dna_apolar'

stored.PDI_Viz_bGTZeroProt = 'bGreaterZeroProtein'
stored.PDI_Viz_qGTZeroProt = 'qGreaterZeroProtein'
stored.PDI_Viz_bGTZeroDNA = 'bGreaterZeroDNA'
stored.PDI_Viz_qGTZeroDNA = 'qGreaterZeroDNA'

stored.PDI_Viz_bGTZeroProt2 = 'bGZeroProtein_m3'
stored.PDI_Viz_qGTZeroProt2 = 'qGZeroProtein_m3'
stored.PDI_Viz_bGTZeroDNA2 = 'bGZeroDNA_g3'
stored.PDI_Viz_qGTZeroDNA2 = 'qGZeroDNA_g3'
stored.PDI_Viz_bGTZeroDNA3 = 'bGZeroDNA_g4'
stored.PDI_Viz_qGTZeroDNA3 = 'qGZeroDNA_g4'

stored.PDI_Viz_prot_obj = ""
stored.PDI_Viz_prot_obj3 = ""
stored.PDI_Viz_dna_obj= ""
stored.PDI_Viz_dna_obj2= ""
stored.PDI_Viz_dna_obj3= ""
stored.PDI_Viz_dna_obj4= ""
stored.PDI_Viz_dna_obj5= ""
stored.PDI_Viz_dna_obj6= ""

stored.PDI_Viz_Threshold = 0.0
stored.PDI_Viz_currentObject = ""

stored.PDI_Viz_Results = {}
stored.PDI_Viz_Stats = {}

#atom removal control
stored.PDI_Viz_RemoveHETATM = 1
stored.PDI_Viz_RemoveH = 1
stored.PDI_Viz_RemoveAltLocs = 1


stored.PDI_Viz_calcmode=0
version = '1.2.3'
name = 'PDIviz'
plugin_name = name+" "+version
def __init__(self):
   self.menuBar.addmenuitem('Plugin', 'command',
                            plugin_name,
                            label = plugin_name,
                            command = lambda s=self : startGUI(s))
def startGUI(app):
  if not stored.PDI_Viz_IsGUIOpen:
    dnavizGUI(app)


about_text = """PDIviz 1.2.1 [2016-04-07]
Analysis and visualization of protein-DNA binding interfaces
Authors: Judemir Ribeiro, Francisco Melo, Andreas Schueller
Contact: aschueller@bio.puc.cl
Website: http://melolab.org/pdiviz

Licensed under the terms of the MIT license. Copyright (c) 2015 Judemir Ribeiro, Francisco Melo, Andreas Schueller

Incorporates the color_b.py script by Robert L. Campbell, Copyright (c) 2004. Used with permission.

PDIviz is a plugin for the PyMOL molecular visualization system that is specifically designed to visualize protein-DNA interfaces in three dimensions. PDIviz detects the interface of protein-DNA complexes by calculating the buried surface area directly in PyMOL. The plugin provides a total of nine distinct 3D visualization modes to highlight interactions with DNA bases and backbone, major and minor groove, and with atoms of different pharmacophoric type.

Usage:
1. Load a protein-DNA complex in PyMOL, e.g. execute the command: fetch 1bl0
2. Open the plugin
3. Click on one of the nine visualization buttons

Original algorithm: Andreas Schueller <aschueller@bio.puc.cl>
Plugin implementation: Judemir Ribeiro <jribeiro@uc.cl>
Project leader: Francisco Melo <fmelo@bio.puc.cl>
"""
Exp1 = """"""

Lend = os.linesep
OSsep = os.sep

#default colors
d_set = [(0.7,0.7,0.7),                                                                   #0: background (RGB)
         (1.0,1.0,1.0),                                                                    #1: DNA no interaction color (RGB)
         (0.6,0.2,0.2),                                                                    #2: DNA backbone interaction (RGB)
         (0.5,0.5,1.0),                                                                   #3: DNA base interaction (RGB)
         (1.0,1.0,0.5),                                                                    #4: DNA both interaction (RGB)
          0.1,                                                                            #5: Transparency (range 0,1)
         (1.0,1.0,1.0, 0.5,0.5,1.0, 0.0,0.0,1.0),                                         #6: grad0 blue (3 RGB colors, low,mid,high)
         (1.0,1.0,1.0, 1.0,0.5,0.5, 1.0,0.0,0.0),                                         #7: grad1 red  (3 RGB colors, low,mid,high)
         (1.0,1.0,1.0, 1.0,1.0,0.5, 1.0,1.0,0.0),                                         #8: grad2 yellow (3 RGB colors, low,mid,high)
         (0.50588,0.50588,0.50588, 1.0,0.0,0.0, 0.0,0.0,1.0, 1.0,1.0,0.0, 1.0,1.0,1.0 )]  #9: rohs visualization colors
                                                                                             #grey for no interactions,
                                                                                             #red H bond acceptor
                                                                                             #blue H bond donor
                                                                                             #yellow thymine methyl
                                                                                             #white carbon with hydrogen
sd_set = [str(item) for item in d_set]


ruby='#993333'
slate='#7F7FFF'
pyellow='#FFFF7F'
iruby='#66CCCC'
islate='#808000'
ipyellow='#000080'



################################################################################
#GUI
class dnavizGUI:


  def GetObjects(self):
    loaded_obj = cmd.get_object_list()
    return [obj for obj in loaded_obj if IsObjProtDNAcomplex(obj) and not IsObjMultiState(obj)]

  def execute(self,result,refocus = True):
    if (result == 'Calculate'):
      cmd.set('suspend_updates', 'on')
      m = '0'
      obj_list = self.selection
      r = Visualize(obj_list[0],m,str(self.set[0]),
            str(self.set[1]),
            str(self.set[2]),
            str(self.set[3]),
            str(self.set[4]),
            str(self.set[5]),
            str(self.set[6]),
            str(self.set[7]),
            str(self.set[8]),
            str(self.set[9]),
            "1")
      cmd.set('suspend_updates', 'off')
      return r
    if (result == 'Color Settings'):
      if not stored.PDI_Viz_IsColorSelOpen:
        color_select(self.parent,self)

    if (result == 'Repaint'):

      self.Repaint()

    if (result == 'Save Image'):
      if not stored.PDI_Viz_IsSaveOpen:
        PNGWrapper(self.parent,self)

    if (result == 'Exit'):
      self.quit()
    if (result == 'Help'):
      if not stored.PDI_Viz_IsHelpOpen:
        self.HelpWindow_p = HelpWindow(self.parent)
    if (result == 'Save Data'):
      self.getpath()

  def Repaint(self):
    self.set[5] = float(self.transp.get())/100.0
    m = '0'
    print("REPAINTING")
    cmd.wizard("message","REPAINTING")

    cmd.set('suspend_updates', 'on')
    obj = self.selection[0]
    Visualize(obj,m,str(self.set[0]),
                 str(self.set[1]),
                  str(self.set[2]),
          str(self.set[3]),
          str(self.set[4]),
          str(self.set[5]),
          str(self.set[6]),
          str(self.set[7]),
          str(self.set[8]),
          str(self.set[9]),
          "0")
    cmd.set('suspend_updates', 'off')
    self.last_mode()
    print("REPAINTING FINSIHED")
    cmd.wizard()
  def quit(self):
    if stored.PDI_Viz_IsColorSelOpen:
      self.c_window.quit()
    stored.PDI_Viz_IsGUIOpen = False
    if stored.PDI_Viz_IsSaveOpen:
      self.PNG_window.quit()
    if stored.PDI_Viz_IsHelpOpen:
      self.HelpWindow_p.quit()
    stored.PDI_Viz_IsSaveOpen = False
    self.dialog.destroy()
  def dummy(self):
    pass


  def __init__(self,app):

    #holds selected object name for the stats page
    self.st= ''
    self.c_window  = ''
    self.PNG_window = ''
    #initializes settings
    self.set = d_set[:]
    self.selection = [""]
    self.parent = app.root
    self.last_mode = self.dummy
    # Create the dialog.
    self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Save Image','Save Data','Exit'),
                                 title = plugin_name,
                                 command = self.execute)
    self.dialog.withdraw()
    self.dialog.protocol('WM_DELETE_WINDOW',self.quit)
    hull = self.dialog.component('hull')
    Pmw.setbusycursorattributes(hull)
    if(os.name == "posix"):
      hull.minsize(440,670)
    elif(os.name == "nt"):
      hull.minsize(402,600)
    hull.resizable(1,1)
    self.notebook = Pmw.NoteBook(self.dialog.interior())
    self.notebook.pack(fill='both',expand=1,padx=10,pady=10)

  # Set up the Main page
    page = self.notebook.add('Main')
    self.group_main = Pmw.Group(page,tag_text='Select Object')
    self.group_main.pack(fill = 'x', expand = 0, padx = 10, pady = 5,anchor='n')
    self.selector = Pmw.OptionMenu(self.group_main.interior(),
                command=self.obj_selector)
    self.group_main.component('hull').bind("<Enter>",self.struct_updater)
    self.selector.pack(fill='x',padx=4,pady=1,side='top')
    group = Pmw.Group(page,tag_text='Toggle Visualizations')
    group.pack(fill = 'x', expand = 0, padx = 10, pady = 5,anchor='n')

    if (os.name == "posix"):
      ###begin X11 block
      self.modebuttons = Pmw.ButtonBox(group.interior()
                         #labelpos = 'nw',
                         #label_text= ':'
                         )
      self.modebuttons.add("Backbone and\n Bases\nSurface",command = self.t1)
      self.modebuttons.add("Major Minor\n Groove\nSurface",command = self.t3)
      self.modebuttons.add("Base\n Interactions\nSurface",command = self.t2)
      self.modebuttons.alignbuttons()
      self.modebuttons2 = Pmw.ButtonBox(group.interior(),
                         labelpos = 'nw',
                         label_text= ' ')
      self.modebuttons2.add("Backbone and\n Bases\nSticks",command = self.t4)
      self.modebuttons2.add("Major Minor\n Groove\nSticks",command = self.t6)
      self.modebuttons2.add("Base\n Interactions\nSticks",command = self.t5)
      self.modebuttons2.alignbuttons()

      self.modebuttons3 = Pmw.ButtonBox(group.interior(),
                         labelpos = 'nw',
                         label_text= ' ')
      self.modebuttons3.add("Backbone and\n Bases\nNucleotides",command = self.t7)
      self.modebuttons3.add("Major Minor\n Groove\nNucleotides",command = self.t9)
      self.modebuttons3.add("Base\n Interactions\nNucleotides", command = self.t8)
      self.modebuttons3.alignbuttons()

      self.modebuttons.grid(column=0,row=0,columnspan=3)
      self.modebuttons2.grid(column=0,row=1,columnspan=3)
      self.modebuttons3.grid(column=0,row=2,columnspan=3)

  #background and transparency buttons
      group3 = Pmw.Group(page,tag_text='Background and transparency')
      group3.pack(fill = 'x', expand = 0, padx = 10, pady = 5)
      self.bgbutton = Pmw.ButtonBox(group3.interior())
      self.bgbutton.add(' Background ',command = self.bg_color,bg = '#B2B2B2',fg = '#000000',height=2,width=10)

      self.transp = Pmw.Counter(group3.interior(),datatype = {'counter' : 'integer'},
              entryfield_value = '10',
              entryfield_validate = {'validator' : 'integer', 'min' : 0,'max':100,'minstrict':1, 'maxstrict':1},
              entryfield_modifiedcommand = self.trans_rt_efield,
              entry_width=3,
              #labelpos = 'nw',
              #label_text = 'Protein surface transparency %:',
              increment=5)
      self.transd = Pmw.Counter(group3.interior(),datatype = {'counter' : 'integer'},
              entryfield_value = '10',
              entryfield_validate = {'validator' : 'integer', 'min' : 0,'max':100,'minstrict':1, 'maxstrict':1},
              entryfield_modifiedcommand = self.trans_rt_efield_d,
              entry_width=3,
              #labelpos = 'nw',
              #label_text = 'DNA surface transparency %:',
              increment=5)
      self.prot_t_lbl = tk.Label(group3.interior(),
                                text = 'Prot Transp: ',
                justify='left',
                anchor='w',
                height=1)
      self.dna_t_lbl = tk.Label(group3.interior(),
                                text = 'DNA Transp: ',
                justify='left',
                anchor='w',
                height=1)
      self.spacer1 = tk.Label(group3.interior(),
                            text = '  ',
                            height = 1)
      self.spacer2 = tk.Label(group3.interior(),
                            text = '           ',
                            height = 1)
      self.spacer3 = tk.Label(group3.interior(),
                            text = '  ',
                            height = 1)
      self.spacer4 = tk.Label(group3.interior(),
                            text = '           ',
                            height = 1)
      self.bgbutton.grid(column=0,row=0,rowspan=3)

      self.spacer1.grid(column=1,row=0,rowspan=1,columnspan=2)
      self.prot_t_lbl.grid(column=3,row=0,rowspan=1,columnspan=2,sticky='W')
      self.spacer2.grid(column=5,row=0,rowspan=1,columnspan=3)
      self.dna_t_lbl.grid(column=8,row=0,rowspan=1,columnspan=2,sticky='W')

      self.spacer3.grid(column=1,row=1,rowspan=1,columnspan=2)

      self.transp.grid(column=3,row=1,rowspan=1,columnspan=2)
      self.spacer4.grid(column=5,row=1,rowspan=1,columnspan=3)
      self.transd.grid(column=8,row=1,rowspan=1)
      ###end X11 block
    elif(os.name == "nt"):
      ###begin Winblows block

      self.modebuttons = Pmw.ButtonBox(group.interior(),
                         #labelpos = 'nw',
                         #label_text= ':'
                         padx=12,
                         pady=8,
                         orient='vertical')
      self.modebuttons.add ("   Backbone and   \n Bases\nSurface",command = self.t1)
      self.modebuttons.add ("   Backbone and   \n Bases\nSticks",command = self.t4)
      self.modebuttons.add ("   Backbone and   \n Bases\nNucleotides",command = self.t7)


      self.modebuttons.alignbuttons()

      self.modebuttons2 = Pmw.ButtonBox(group.interior(),
                         #labelpos = 'nw',
                         #label_text= ':'
                         padx=12,
                         pady=8,
                         orient='vertical')
      self.modebuttons2.add("    Major Minor   \n Groove\nSurface",command = self.t3)
      self.modebuttons2.add("    Major Minor   \n Groove\nSticks",command = self.t6)
      self.modebuttons2.add("    Major Minor   \n Groove\nNucleotides",command = self.t9)

      self.modebuttons2.alignbuttons()

      self.modebuttons3 = Pmw.ButtonBox(group.interior(),
                         #labelpos = 'nw',
                         #label_text= ':'
                         padx=12,
                         pady=8,
                         orient='vertical')

      self.modebuttons3.add("Base\n    Interactions    \nSurface",command = self.t2)
      self.modebuttons3.add("Base\n    Interactions    \nSticks",command = self.t5)
      self.modebuttons3.add("Base\n    Interactions    \nNucleotides", command = self.t8)

      self.modebuttons3.alignbuttons()

      self.modebuttons.grid(column=0,row=0,rowspan=3)
      self.modebuttons2.grid(column=1,row=0,rowspan=3)
      self.modebuttons3.grid(column=2,row=0,rowspan=3)

  #background and transparency buttons
      group3 = Pmw.Group(page,tag_text='Background and transparency')
      group3.pack(fill = 'x', expand = 0, padx = 10, pady = 5)
      self.bgbutton = Pmw.ButtonBox(group3.interior())
      self.bgbutton.add('Background',command = self.bg_color,bg = '#B2B2B2',fg = '#000000',height=3,width=14)

      self.transp = Pmw.Counter(group3.interior(),datatype = {'counter' : 'integer'},
              entryfield_value = '10',
              entryfield_validate = {'validator' : 'integer', 'min' : 0,'max':100,'minstrict':1, 'maxstrict':1},
              entryfield_modifiedcommand = self.trans_rt_efield,
              entry_width=3,
              #labelpos = 'nw',
              #label_text = 'Protein surface transparency %:',
              increment=5)
      self.transd = Pmw.Counter(group3.interior(),datatype = {'counter' : 'integer'},
              entryfield_value = '10',
              entryfield_validate = {'validator' : 'integer', 'min' : 0,'max':100,'minstrict':1, 'maxstrict':1},
              entryfield_modifiedcommand = self.trans_rt_efield_d,
              entry_width=3,
              #labelpos = 'nw',
              #label_text = 'DNA surface transparency %:',
              increment=5)
      self.prot_t_lbl = tk.Label(group3.interior(),
                                text = 'Prot Transp: ',
                justify='left',
                anchor='w',
                height=1)
      self.dna_t_lbl = tk.Label(group3.interior(),
                                text = 'DNA Transp: ',
                justify='left',
                anchor='w',
                height=1)

      sp13 = 2
      sp24 = 7
      self.spacer1 = tk.Label(group3.interior(),
                            text = ' ',
                            height = 1,
                            width = sp13)
      self.spacer2 = tk.Label(group3.interior(),
                            text = ' ',
                            height = 1,
                            width = sp24)
      self.spacer3 = tk.Label(group3.interior(),
                            text = ' ',
                            height = 1,
                            width = sp13)
      self.spacer4 = tk.Label(group3.interior(),
                            text = ' ',
                            height = 1,
                            width = sp24)
      self.spacer5 = tk.Label(group3.interior(),
                            text = ' ',
                            height = 1)
      self.spacer5.grid(column=0,row=0,rowspan=3)
      self.bgbutton.grid(column=1,row=0,rowspan=3)

      self.spacer1.grid(column=2,row=0,rowspan=1,columnspan=2)
      self.prot_t_lbl.grid(column=4,row=0,rowspan=1,columnspan=2,sticky='W')
      self.spacer2.grid(column=6,row=0,rowspan=1,columnspan=3)
      self.dna_t_lbl.grid(column=9,row=0,rowspan=1,columnspan=2,sticky='W')

      self.spacer3.grid(column=2,row=1,rowspan=1,columnspan=2)

      self.transp.grid(column=4,row=1,rowspan=1,columnspan=2)
      self.spacer4.grid(column=6,row=1,rowspan=1,columnspan=3)
      self.transd.grid(column=9,row=1,rowspan=1)
      ###end Windows block
    spopt = Pmw.Group(page,tag_text='Special options')
    spopt.pack(fill = 'x', expand = 0, padx = 10, pady = 5,anchor='n')
    explanation = tkinter.Label(spopt.interior(),
                                text = "BSA values greater than this are considered as interaction.\nPRESS ENTER TO REPAINT.",
                justify='left',
                anchor='w',
                height=3
                               )
    explanation.pack(fill='x',padx=4,pady=1,expand=0)



    self.cutoff = Pmw.EntryField(spopt.interior(),
                   validate = { 'validator' : 'real',
                                  'min' : 0,
                                  'minstrict': 0},
                 modifiedcommand = self.set_cutoff,
                 command = self.Repaint,
                 value = '0.0'
                 )
    self.cutoff.pack(fill='x',padx=4,pady=1,side='top')
  #set up stats page
    page = self.notebook.add('Statistics')
    self.group1 = Pmw.Group(page,tag_text = 'Stats:')
    self.group1.pack(fill = 'both', expand = 0, padx = 4, pady = 5)
    self.table = SimpleTable(self.group1.interior(),28,2)
    self.table.pack()
    spacer = tkinter.Label(self.group1.interior(),
                                text = '',
                justify='left',
                anchor='w',
                height=1
                               )
    spacer.pack()
    tkinter.Button(self.group1.interior(),text='Copy to clipboard', command=self.copystat,height=3).pack()

    #first column labels
    self.labels = ['Molecular Object:','Complex SASA','Free protein SASA','Free DNA SASA','',
    'Protein backbone SASA','Protein side chain SASA','DNA backbone SASA','DNA bases SASA',
    'DNA major groove SASA','DNA minor groove SASA','','Buried protein surface','Buried protein backbone surface','Buried protein side chain surface','',
    'Buried DNA surface','Buried DNA backbone surface','Buried DNA bases surface','Buried DNA major groove surface','Buried DNA minor groove surface','',
    'Buried Protein polar surface','Buried Protein apolar surface','Buried DNA polar surface','Buried DNA apolar surface','',
    'Interface Area']
    for i in range(len(self.labels)):
      self.table.set(i,0,self.labels[i])
    self.table.set(0,1,'Area [A^2]')

  #set up about page
    page = self.notebook.add('About')
    group = Pmw.Group(page, tag_text='About '+plugin_name)
    group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
    text = about_text
    interior_frame=tk.Frame(group.interior())
    bar=tk.Scrollbar(interior_frame)
    text_holder=tk.Text(interior_frame,yscrollcommand=bar.set,background="#ddddff",font="Times 10",wrap='word')
    bar.config(command=text_holder.yview)
    text_holder.insert('end',text)
    text_holder.config(state='disabled')
    text_holder.pack(side='left',expand="yes",fill="both")

    bar.pack(side='left',expand="yes",fill="y")
    interior_frame.pack(expand="yes",fill="both")

  #initialize selector
    self.obj_list = self.GetObjects()
    if len(self.obj_list) > 0:
      self.selector.setitems(self.obj_list)
      self.selector.setvalue(self.obj_list[0])
      self.obj_selector(self.obj_list[0])
  #show dialog and set opened window flag

    self.dialog.show()
    stored.PDI_Viz_IsGUIOpen = True

#support functions
  def struct_updater(self,event):

    obj_list = self.GetObjects()
    if obj_list != self.obj_list:
      self.obj_list = obj_list
      self.selector.destroy()
      self.selector = Pmw.OptionMenu(self.group_main.interior(),
              command=self.obj_selector)
      self.selector.pack(fill='x',padx=4,pady=1,side='top')
      self.selector.setitems(self.obj_list)
      if len(self.obj_list) > 0:
        self.selector.setvalue(self.obj_list[0])
        self.obj_selector(self.obj_list[0])

  def obj_selector(self,obj_name):
    obj_list = cmd.get_object_list()
    for obj in obj_list:
      if(obj_name == obj):
        cmd.enable(obj)
      else:
        cmd.disable(obj)
    self.selection = [obj_name]
    self.statselector(obj_name)
    cmd.origin(obj_name)
    cmd.center(obj_name)
    self.sinkmodebuttons(-1,-1)
    #self.last_mode(recalc=False)
    #self.Repaint()
  def byrescallback(self,tag,state):
    if state:
      stored.bb_ba_byres = True
    else:
      stored.bb_ba_byres = False
    self.Repaint()

  def sinkmodebuttons(self,x,y):
    modebuttons = [self.modebuttons,self.modebuttons2,self.modebuttons3]
    if(os.name == "posix"):
      for row in modebuttons:
        for i in range(row.numbuttons()):
          row.button(i).config(relief='raised',state='normal')
      if(x >= 0 and y >= 0):
        modebuttons[y].button(x).config(relief='sunken',state='active')
    if(os.name == "nt"):
      for col in modebuttons:
        for i in range(col.numbuttons()):
          col.button(i).config(relief='raised',state='normal')
      if(x >= 0 and y >= 0):
        modebuttons[x].button(y).config(relief='sunken',state='active')


  def t1(self,recalc=True):
    self.sinkmodebuttons(0,0)
    cur = ["_dna_obj_m3","_prot_obj_m2"]
    oth = ["_dna_int_m2","_dna_obj_g3","_prot_int_m2","_prot_obj_g3","_dna_obj_m2","_dna_int_m3","_dna_obj_g4"]
    self.last_mode = self.t1
    self.toggle1(0,cur,oth,recalc)
  def t2(self,recalc=True):
    self.sinkmodebuttons(2,0)
    cur = ["_dna_int_m2","_prot_int_m2"]
    oth = ["_dna_obj_m3","_dna_obj_g3","_prot_obj_g3","_prot_obj_m2","_dna_obj_m2","_dna_int_m3","_dna_obj_g4"]
    self.last_mode = self.t2
    self.toggle1(0,cur,oth,recalc)
  def t3(self,recalc=True):
    self.sinkmodebuttons(1,0)
    cur = ["_dna_obj_g3","_prot_obj_g3"]
    oth = ["_dna_obj_m3","_dna_int_m2","_prot_obj_m2","_prot_int_m2","_dna_obj_m2","_dna_int_m3","_dna_obj_g4"]
    self.last_mode = self.t3
    self.toggle1(0,cur,oth,recalc)
  def t4(self,recalc=True):
    self.sinkmodebuttons(0,1)
    cur = ["_dna_obj_m3","_prot_obj_m2"]
    oth = ["_dna_int_m2","_dna_obj_g3","_prot_int_m2","_prot_obj_g3","_dna_obj_m2","_dna_int_m3","_dna_obj_g4"]
    self.last_mode = self.t4
    self.toggle1(1,cur,oth,recalc)
  def t5(self,recalc=True):
    self.sinkmodebuttons(2,1)
    cur = ["_dna_int_m2","_prot_int_m2"]
    oth = ["_dna_obj_m3","_dna_obj_g3","_prot_obj_g3","_prot_obj_m2","_dna_obj_m2","_dna_int_m3","_dna_obj_g4"]
    self.last_mode = self.t5
    self.toggle1(1,cur,oth,recalc)
  def t6(self,recalc=True):
    self.sinkmodebuttons(1,1)
    cur = ["_dna_obj_g3","_prot_obj_g3"]
    oth = ["_dna_obj_m3","_dna_int_m2","_prot_obj_m2","_prot_int_m2","_dna_obj_m2","_dna_int_m3","_dna_obj_g4"]
    self.last_mode = self.t6
    self.toggle1(1,cur,oth,recalc)
  def t7(self,recalc=True):
    self.sinkmodebuttons(0,2)
    cur = ["_dna_obj_m2","_prot_obj_m2"]
    oth = ["_dna_int_m2","_dna_obj_g3","_prot_int_m2","_prot_obj_g3","dna_obj_m3","_dna_int_m3","_dna_obj_g4"]
    self.last_mode = self.t7
    self.toggle1(1,cur,oth,recalc)
  def t8(self,recalc=True):
    self.sinkmodebuttons(2,2)
    cur = ["_dna_int_m3","_prot_int_m2"]
    oth = ["_dna_int_m2","_dna_obj_g3","_prot_obj_m2","_prot_obj_g3","dna_obj_m3","_dna_obj_m2","_dna_obj_g4"]
    self.last_mode = self.t8
    self.toggle1(1,cur,oth,recalc)
  def t9(self,recalc=True):
    self.sinkmodebuttons(1,2)
    cur = ["_dna_obj_g4","_prot_obj_g3"]
    oth = ["_dna_int_m2","_dna_obj_g3","_prot_int_m2","_prot_obj_m2","dna_obj_m3","_dna_int_m3","_dna_obj_m2"]
    self.last_mode = self.t9
    self.toggle1(1,cur,oth,recalc)

  def toggle1(self,m,cur,oth,recalc=True):
    obj = cmd.get_object_list()
    sele = self.selection
    sobj = []
    currentobj = []
    if sele[0] == "":
      self.struct_updater(0)
      if self.selection[0] != "":
        self.obj_selector(self.selection[0])
        sele = self.selection
      else:
        self.sinkmodebuttons(-1,-1)
        return


    if recalc:
      for item in [sele[0]+cur[0],sele[0]+cur[1]]:
        if not item in obj:
          print("Please wait...")
          cmd.wizard("message","Please wait...")
          if(not self.execute('Calculate')):
            cmd.wizard()
            self.sinkmodebuttons(-1,-1)
            return
          else:
            print("CALCULATION FINISHED...")
            self.statselector(sele[0])
          #self.toggle1(m,cur,oth)
            cmd.wizard()
            break
          #return
    cmd.wizard()
    obj = cmd.get_object_list()
    for item in obj:
      for s in sele:
        cmd.disable(s)
        if s in item:
          sobj.append(item)
    for item in sobj:
      for c in cur:
        if c in item:
          cmd.enable(item)
          while(not item in cmd.get_names(enabled_only=1)):
            cmd.enable(item)
          currentobj.append(item)
    for item in sobj:
      for o in oth:
        if o in item:
          cmd.disable(item)
    if m == 0: #surface
      for item in currentobj:
        if cur[0] in item:
          cmd.show('surface',item)
          cmd.hide('cartoon',item)
        if cur[1] in item:
          cmd.hide('surface',item)
    else:
      for item in currentobj:
        if cur[1] in item:
          cmd.show('surface',item)
        if cur[0] in item:
          cmd.hide('surface',item)
          cmd.show('cartoon',item)
    cmd.wizard()
  def set_cutoff(self):
    try:
      val = float(self.cutoff.getvalue())
      stored.PDI_Viz_dSASAcutoff = val
    except:
      stored.PDI_Viz_dSASAcutoff = 0;

  def trans_rt_efield(self):
    otype = ["_prot_obj_m2","_prot_obj_g3","_prot_int_m2"]
    sele = self.selection
    objects = cmd.get_object_list()
    objects = [obj for obj in objects if otype[0] in obj or otype[1] in obj or otype[2] in obj]
    sobjects = []
    for obj in objects:
      for s in sele:
        if s in obj:
          sobjects.append(obj)
    text = self.transp.component('entryfield').getvalue()
    try:
      v = int(text)
    except:
      return
    for obj in sobjects:
      cmd.set("transparency",str(v/100.0),"object "+obj)

  def trans_rt_efield_d(self):
    otype = ["_dna_obj_m3","_dna_int_m2","_dna_obj_g3"]

    sele = self.selection
    objects = cmd.get_object_list()
    objects = [obj for obj in objects if otype[0] in obj or otype[1] in obj or otype[2] in obj]
    sobjects = []
    for obj in objects:
      for s in sele:
        if s in obj:
          sobjects.append(obj)
    text = self.transd.component('entryfield').getvalue()
    try:
      v = int(text)
    except:
      return
    for obj in sobjects:
      cmd.set("transparency",str(v/100.0),"object "+obj)

  def statselector(self,obj_name):
    try:
      stats = stored.PDI_Viz_Stats[obj_name]
      self.group1.component('tag').config(text = "Stats from "+obj_name+":")
    except:
      for i in range(1,28):
        self.table.set(i,1,'')
      return
    self.st = obj_name
    self.table.set(1,1,stats['total'])
    self.table.set(2,1,stats['protein'])
    self.table.set(3,1,stats['dna'])

    self.table.set(5,1,stats['prot_bb'])
    self.table.set(6,1,stats['prot_sc'])
    self.table.set(7,1,stats['dna_bb'])
    self.table.set(8,1,stats['dna_ba'])
    self.table.set(9,1,stats['dna_ma'])
    self.table.set(10,1,stats['dna_mi'])

    self.table.set(12,1,stats['prot-total'])
    self.table.set(13,1,stats['prot-total_bb'])
    self.table.set(14,1,stats['prot-total_sc'])

    self.table.set(16,1,stats['dna-total'])
    self.table.set(17,1,stats['dna-total_bb'])
    self.table.set(18,1,stats['dna-total_ba'])
    self.table.set(19,1,stats['dna-total_ma'])
    self.table.set(20,1,stats['dna-total_mi'])
    
    self.table.set(22,1,stats['prot-total_prot_polar'])
    self.table.set(23,1,stats['prot-total_prot_apolar'])
    self.table.set(24,1,stats['dna-total_dna_polar'])
    self.table.set(25,1,stats['dna-total_dna_apolar'])
    
    self.table.set(27,1,stats['interface'])

  def copystat(self):
    if self.st in stored.PDI_Viz_Stats:
      stats = stored.PDI_Viz_Stats[self.st]
      labels = self.labels
      res = labels[0]+'\t'+'Area [A^2]' + Lend
      res+=labels[1]+'\t'+ stats['total'] + Lend
      res+=labels[2]+'\t'+ stats['protein'] + Lend
      res+=labels[3]+'\t'+ stats['dna'] + Lend + Lend
      res+=labels[5]+'\t'+ stats['prot_bb'] + Lend
      res+=labels[6]+'\t'+ stats['prot_sc'] + Lend
      res+=labels[7]+'\t'+ stats['dna_bb'] + Lend
      res+=labels[8]+'\t'+ stats['dna_ba'] + Lend
      res+=labels[9]+'\t'+ stats['dna_ma'] + Lend
      res+=labels[10]+'\t'+ stats['dna_mi'] + Lend + Lend
      res+=labels[12]+'\t'+ stats['prot-total'] + Lend
      res+=labels[13]+'\t'+ stats['prot-total_bb'] + Lend
      res+=labels[14]+'\t'+ stats['prot-total_sc'] +Lend+Lend
      res+=labels[16]+'\t'+ stats['dna-total'] + Lend
      res+=labels[17]+'\t'+ stats['dna-total_bb'] + Lend
      res+=labels[18]+'\t'+ stats['dna-total_ba'] + Lend
      res+=labels[19]+'\t'+ stats['dna-total_ma'] + Lend
      res+=labels[20]+'\t'+ stats['dna-total_mi'] + Lend + Lend 
      res+=labels[22]+'\t'+ stats['prot-total_prot_polar'] + Lend 
      res+=labels[23]+'\t'+ stats['prot-total_prot_apolar'] + Lend 
      res+=labels[24]+'\t'+ stats['dna-total_dna_polar'] + Lend 
      res+=labels[25]+'\t'+ stats['dna-total_dna_apolar'] + Lend + Lend
      res+=labels[27]+'\t'+ stats['interface'] + Lend

      r = tk.Tk()
      r.withdraw()
      r.clipboard_clear()
      r.clipboard_append(res)
      r.destroy()

  def getpath(self):
    path = tkinter.filedialog.askdirectory(title='Select path',mustexist=False)
    #print type(path),"|"+path+"|"
    if (type(path) == type(str()) or type(path) == type(str())) and path !="" :
      path = os.path.normpath(path)
      path += OSsep
      if self.selection[0]:
        PrintOutputFiles(self.selection[0],path)

  def modechange(self,tag):
    if stored.PDI_Viz_IsColorSelOpen:
      cwin = self.c_window
      if tag == 'Mode 1':
        cwin.grad0.component('label').config(text = 'DNA contact area:',justify = 'left')
        cwin.grad1.component('label').config(text = 'Disabled in mode 1:',justify = 'left')
        cwin.grad2.component('label').config(text = 'Disabled in mode 1:',justify = 'left')
        g1b = cwin.grad1.button
        g2b = cwin.grad2.button
        for i in range(3):
          g1b(i).config(state='disabled')
          g2b(i).config(state='disabled')

      if tag == 'Mode 2':
        cwin.grad0.component('label').config(text = 'Gradient 1:\nDNA base contact:',justify = 'left')
        cwin.grad1.component('label').config(text = 'Gradient 2:\nDNA backbone contact:',justify = 'left')
        cwin.grad2.component('label').config(text = 'Gradient 3:\nDNA base and backbone contact:',justify = 'left')
        g1b = cwin.grad1.button
        g2b = cwin.grad2.button
        for i in range(3):
          g1b(i).config(state='normal')
          g2b(i).config(state='normal')

  def bg_color(self):
    b = self.bgbutton.button(0)
    bg = b.config()['background'][4]
    res = tkinter.colorchooser.askcolor(bg,title= 'Background color')
    if(res[0] != None):
      r = res[0]
      s = [r[0]/255.0,r[1]/255.0,r[2]/255.0]
      self.set[0] = tuple(s)
      ri = '#%02x%02x%02x' % (255-r[0],255-r[1],255-r[2])
      b.config(bg=res[1],fg=ri)
      cmd.set("bg_rgb",str(s))

#stat table class
class SimpleTable(tk.Frame):
  def __init__(self, parent, rows=10, columns=2):
        # use black background so it "peeks through" to
        # form grid lines
    tk.Frame.__init__(self, parent, background="black")
    self._widgets = []
    for row in range(rows):
      current_row = []
      for column in range(columns):
        label = tk.Label(self, text='',borderwidth=0, width=10)
        label.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
        current_row.append(label)
      self._widgets.append(current_row)

  def set(self, row, column, value):
    widget = self._widgets[row][column]
    if (column == 0):
      widget.configure(text=value,width=30,justify='left',anchor='w')
    else:
      if(row!=0):
        widget.configure(text=value,width=15,justify='right',anchor='e')
      else:
        widget.configure(text=value,width=15)

#color selection window
class color_select:
  def execute(self,result,refocus = True):

    if (result == 'Repaint'):
      self.mw.Repaint()
    if (result == 'Default Settings'):
      self.mw.set = d_set[:]
      self.dnabutton.button(0).config(bg='#FFFFFF',fg='#000000')
      self.dnabutton.button(1).config(bg=slate,fg=islate)
      self.dnabutton.button(2).config(bg=ruby,fg=iruby)
      self.dnabutton.button(3).config(bg=pyellow,fg=ipyellow)
      self.dnabutton2.button(0).config(bg='#818181',fg='#7E7E7E')
      self.dnabutton2.button(1).config(bg='#FF0000',fg='#00FFFF')
      self.dnabutton2.button(2).config(bg='#0000FF',fg='#0000FF')
      self.dnabutton2.button(3).config(bg='#FFFF00',fg='#0000FF')
      self.dnabutton2.button(4).config(bg='#FFFFFF',fg='#000000')
      self.grad0.button(0).config(bg='#FFFFFF',fg='#000000')
      self.grad0.button(1).config(bg='#8080FF',fg='#7F7F00')
      self.grad0.button(2).config(bg='#0000FF',fg='#FFFF00')
      self.grad1.button(0).config(bg='#FFFFFF',fg='#000000')
      self.grad1.button(1).config(bg='#FF8080',fg='#007F7F')
      self.grad1.button(2).config(bg='#FF0000',fg='#00FFFF')
      self.grad2.button(0).config(bg='#FFFFFF',fg='#000000')
      self.grad2.button(1).config(bg='#FFFF80',fg='#00007F')
      self.grad2.button(2).config(bg='#FFFF00',fg='#0000FF')
      self.mw.transp.setentry(10)
      self.mw.bgbutton.button(0).config(bg = '#808080',fg = '#000000')
    if (result == 'Close'):
      self.quit()

  def quit(self):
    stored.PDI_Viz_IsColorSelOpen = False
    self.dialog.destroy()

  def __init__(self,app,main_w):
    #pointers to main window

    self.parent = app
    self.mw = main_w
    self.mw.set = d_set[:]

    self.mw.c_window = self
    # Create the dialog.
    self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Repaint','Default Settings' ,'Close'),
                                 title = plugin_name + ' - Color Settings',
                                 command = self.execute)
    self.dialog.withdraw()
    self.dialog.protocol('WM_DELETE_WINDOW',self.quit)
    hull = self.dialog.component('hull')
    Pmw.setbusycursorattributes(hull)
    #dialog window size
    hull.minsize(740,580)
    hull.resizable(0,0)
    self.notebook = Pmw.NoteBook(self.dialog.interior())
    self.notebook.pack(fill='both',expand=1,padx=10,pady=10)
  #create Backbone and Bases page
    page = self.notebook.add('Backbone and\nBases')
    #page = self.dialog.interior()


    group1 = Pmw.Group(page,tag_text='DNA nucleotide colors')
    group1.pack(fill = 'x', expand = 0, padx = 4, pady = 5)
    #dna color buttons
    self.dnabutton = Pmw.ButtonBox(group1.interior(),labelpos= 'nw')
    self.dnabutton.pack(fill= 'x', expand = 0)
    self.dnabutton.add('No\nInteraction',command = self.dna_normal,bg='#FFFFFF',fg='#000000',height=3,width=12)
    self.dnabutton.add('Base\nInteraction',command = self.dna_base,bg=slate,fg=islate,height=3,width=12)
    self.dnabutton.add('Backbone\nInteraction',command = self.dna_bb,bg=ruby,fg=iruby,height=3,width=12)
    self.dnabutton.add('Base and\nBackbone\nInteraction',command = self.dna_both,bg=pyellow,fg=ipyellow,height=3,width=12)

    #protein gradient buttons
    group2 = Pmw.Group(page,tag_text='Protein surface gradients')
    group2.pack(fill = 'x', expand = 0, padx = 4, pady = 5)

    self.grad0 = Pmw.ButtonBox(group2.interior(),labelpos= 'nw',label_text = 'DNA contact area:')
    self.grad0.pack(fill = 'x', expand = 0)
    self.grad0.add('Low',command = self.grad0_low,bg='#FFFFFF',fg='#000000',height=3,width=10)
    self.grad0.add('Mid',command = self.grad0_mid,bg='#8080FF',fg='#7F7F00',height=3,width=10)
    self.grad0.add('High',command = self.grad0_high,bg='#0000FF',fg='#FFFF00',height=3,width=10)

    self.grad1 = Pmw.ButtonBox(group2.interior(),labelpos= 'nw',label_text = 'Disabled in mode 1:')
    self.grad1.pack(fill = 'x', expand = 0)
    self.grad1.add('Low',command = self.grad1_low,state='disabled',bg='#FFFFFF',fg='#000000',height=3,width=10)
    self.grad1.add('Mid',command = self.grad1_mid,state='disabled',bg='#FF8080',fg='#007F7F',height=3,width=10)
    self.grad1.add('High',command = self.grad1_high,state='disabled',bg='#FF0000',fg='#00FFFF',height=3,width=10)


    self.grad2 = Pmw.ButtonBox(group2.interior(),labelpos= 'nw',label_text = 'Disabled in mode 1:')
    self.grad2.pack(fill = 'x', expand = 0)
    self.grad2.add('Low',command = self.grad2_low,state='disabled',bg='#FFFFFF',fg='#000000',height=3,width=10)
    self.grad2.add('Mid',command = self.grad2_mid,state='disabled',bg='#FFFF80',fg='#00007F',height=3,width=10)
    self.grad2.add('High',command = self.grad2_high,state='disabled',bg='#FFFF00',fg='#0000FF',height=3,width=10)
  #create Bases page:
    page = self.notebook.add('Base interactions')
    group3 = Pmw.Group(page,tag_text='DNA chemical bases colors')
    group3.pack(fill = 'x', expand = 0, padx = 4, pady = 5)
    self.dnabutton2 = Pmw.ButtonBox(group3.interior(),labelpos= 'nw')
    self.dnabutton2.pack(fill= 'x', expand = 0)
    self.dnabutton2.add('No\nInteraction',command = self.dna_noint,bg='#818181',fg='#7E7E7E',height=3,width=12)
    self.dnabutton2.add('H bond\nAcceptor',command = self.dna_acceptor,bg='#FF0000',fg='#00FFFF',height=3,width=12)
    self.dnabutton2.add('H bond\nDonor',command = self.dna_donor,bg='#0000FF',fg='#FFFF00',height=3,width=12)
    self.dnabutton2.add('Thymine\nMethyl',command = self.dna_thymine,bg='#FFFF00',fg='#0000FF',height=3,width=12)
    self.dnabutton2.add('Neutral\nCarbon',command = self.dna_neutral,bg='#FFFFFF',fg='#000000',height=3,width=12)

    #Show dialog and set opened dialog flag
    stored.PDI_Viz_IsColorSelOpen = True
    self.mw.modechange("Mode 2")
    self.dialog.show()

  def dna_color(self,button,var,name):

    b = self.dnabutton.button(button)
    bg = b.config()['background'][4]
    res = tkinter.colorchooser.askcolor(bg,title= name)
    if(res[0] != None):
      r = res[0]
      self.mw.set[var] = (r[0]/255.0,r[1]/255.0,r[2]/255.0)
      ri = '#%02x%02x%02x' % (255-r[0],255-r[1],255-r[2])
      b.config(bg=res[1],fg=ri)

  def dna_color2(self,button,var,name):

    b = self.dnabutton2.button(button)
    bg = b.config()['background'][4]
    res = tkinter.colorchooser.askcolor(bg,title= name)
    if(res[0] != None):
      r = res[0]
      copy = list(self.mw.set[var])
      copy[  button*3] = r[0]/255.0
      copy[1+button*3] = r[1]/255.0
      copy[2+button*3] = r[2]/255.0
      self.mw.set[var] = tuple(copy)
      ri = '#%02x%02x%02x' % (255-r[0],255-r[1],255-r[2])
      b.config(bg=res[1],fg=ri)

  def prot_color(self,grad,level,var,name):
    bbox_d = {0:self.grad0,1:self.grad1,2:self.grad2}
    b = bbox_d[grad].button(level)
    bg = b.config()['background'][4]
    res = tkinter.colorchooser.askcolor(bg,title= name)
    if(res[0] != None):
      r = res[0]
      t = list(self.mw.set[var])
      t[0+3*level:3+3*level] = (r[0]/255.0,r[1]/255.0,r[2]/255.0)
      self.mw.set[var] = tuple(t)
      ri = '#%02x%02x%02x' % (255-r[0],255-r[1],255-r[2])
      b.config(bg=res[1],fg=ri)

  def dna_noint(self):
    self.dna_color2(0,9,'No Interaction')
  def dna_acceptor(self):
    self.dna_color2(1,9,'H bond Acceptor')
  def dna_donor(self):
    self.dna_color2(2,9,'H bond Donor')
  def dna_thymine(self):
    self.dna_color2(3,9,'Thymine Methyl')
  def dna_neutral(self):
    self.dna_color2(4,9,'Neutral Carbon')

  def dna_normal(self):
    self.dna_color(0,1,'No DNA interaction')
  def dna_bb(self):
    self.dna_color(2,2,'Backbone interaction')
  def dna_base(self):
    self.dna_color(1,3,'Base interaction')
  def dna_both(self):
    self.dna_color(3,4,'Backbone and base interaction')

  def grad0_low(self):
    self.prot_color(0,0,6,'Gradient 1 Low')
  def grad0_mid(self):
    self.prot_color(0,1,6,'Gradient 1 Mid')
  def grad0_high(self):
    self.prot_color(0,2,6,'Gradient 1 High')

  def grad1_low(self):
    self.prot_color(1,0,7,'Gradient 2 Low')
  def grad1_mid(self):
    self.prot_color(1,1,7,'Gradient 2 Mid')
  def grad1_high(self):
    self.prot_color(1,2,7,'Gradient 2 High')

  def grad2_low(self):
    self.prot_color(2,0,8,'Gradient 3 Low')
  def grad2_mid(self):
    self.prot_color(2,1,8,'Gradient 3 Mid')
  def grad2_high(self):
    self.prot_color(2,2,8,'Gradient 3 High')

class PNGWrapper:
  def __init__(self,app,main_w):
    self.res = (2400,2400)
    self.dpi = 300.0
    self.ray = 1
    self.mode = 1
    self.saving = 0
    self.parent = app
    self.mw = main_w

    self.mw.PNG_window = self
    self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Save...', 'Close'),
                                 title = plugin_name + ' - Save Image...',
                                 command = self.execute)
    self.dialog.withdraw()
    self.dialog.protocol('WM_DELETE_WINDOW',self.quit)
    hull = self.dialog.component('hull')
    Pmw.setbusycursorattributes(hull)
    hull.minsize(512,384)
    hull.resizable(0,0)
    page = self.dialog.interior()
    self.reslist = ["640x480","1024x768","1600x1200","2048x1536","854x480",
                "1280x720", "1920x1080", "2560x1440","3840x2160","768x768",
                "1024x1024","2048x2048","2400x2400","4800x4800"]
    self.resmap = {"640x480":   (640, 480),
                 "1024x768": (1024, 768),
               "1600x1200":(1600,1200),
               "2048x1536":(2048,1536),
               "854x480":   (854, 480),
               "1280x720": (1280, 720),
               "1920x1080":(1920,1080),
               "2560x1440":(2560,1440),
               "3840x2160":(3840,2160),
               "768x768":   (768, 768),
               "1024x1024":(1024,1024),
               "2048x2048":(2048,2048),
               "2400x2400":(2400,2400),
               "4800x4800":(4800,4800)}

    self.dpilist = ["96", "150", "200", "300", "600"]
    self.dpimap = {"96": 96.0, "150": 150.0, "200": 200.0, "300": 300.0, "600": 600.0}

    self.modes = ["Normal Ray Tracing",
                  "Ray Tracing + Black outlines",
                  "Cel Shaded Ray Tracing"]

  #resolution selection
    self.sel_group = Pmw.Group(page,tag_text='Output:')
    self.sel_group.pack(fill = 'x', expand = 0, padx = 4, pady = 5)

    self.selector = Pmw.OptionMenu(self.sel_group.interior(),labelpos = 'nw',
                label_text = 'Select Resolution:',
                command=self.resolution)
    self.selector.pack(fill = 'x', expand = 0, padx = 4, pady = 5)
    self.selector.setitems(self.reslist)
    self.selector.setvalue(self.reslist[12])
  #dpi selection
    self.selector2 = Pmw.OptionMenu(self.sel_group.interior(),labelpos = 'nw',
                label_text = 'Select DPI:',
                command=self.dpiset )
    self.selector2.pack(fill = 'x', expand = 0, padx = 4, pady = 5)
    self.selector2.setitems(self.dpilist)
    self.selector2.setvalue(self.dpilist[3])

  #ray tracing options
    group1 = Pmw.Group(page,tag_text='Image Options')
    group1.pack(fill = 'x', expand = 0, padx = 4, pady = 5)
    #ray tracing enable/disable

    self.RayTrace = Pmw.RadioSelect(group1.interior(),
                     buttontype = 'checkbutton',
                     orient = 'vertical',
                           command = self.checkbox)
    self.RayTrace.pack(fill = 'x')
    for option in ["Enable Ray Tracing","Ray Tracing Self Shadowing"]:
      self.RayTrace.add(option)

    #ray mode
    self.raymode = Pmw.RadioSelect(group1.interior(),
                     buttontype= 'radiobutton',
                       labelpos= 'w',
                       orient = 'vertical',
                       command = self.raymode
                     #label_text = "Ray tracing\n  mode"
                     )
    self.raymode.pack(fill= 'x')
    for mode in self.modes:
      self.raymode.add(mode)
    self.raymode.invoke(self.modes[0])
    self.RayTrace.invoke("Enable Ray Tracing")
    self.RayTrace.invoke("Ray Tracing Self Shadowing")
  ###create dialog
    stored.PDI_Viz_IsSaveOpen = True
    self.dialog.show()

  #support functions
  def execute(self,result):
    if result == 'Close':
      self.quit()
    if result == 'Save...':
      self.Save()

  def quit(self):
    stored.PDI_Viz_IsSaveOpen = False
    self.dialog.destroy()

  def resolution(self,string):
    self.res = self.resmap[string]

  def Save(self):
    if self.saving == 1:
      return
    self.saving = 1
    self.dialog.component('buttonbox').button(0).config(state='disable')
    filename = tkinter.filedialog.asksaveasfilename(defaultextension='.png',
                            filetypes=[('PNG image','.png')],
                            initialfile="image.png",
                            title="Save Image...",
                            parent=self.parent)


    if filename:
      print("Saving Image with filename: "+filename)
      #print "Resolution: ",self.res[0],"x",self.res[1]
      #print "DPI: ",self.dpi
      #print "Ray: ",self.ray
      cmd.wizard("message","Saving image...")

      if self.ray:
        cmd.png(filename,self.res[0],self.res[1],self.dpi,1)
      else:
        cmd.png(filename,self.res[0],self.res[1],self.dpi,0)
        cmd.png(filename,dpi=self.dpi)
      cmd.wizard()
    self.dialog.component('buttonbox').button(0).config(state='active')
    self.saving = 0

  def dpiset(self,dpi):
    self.dpi = self.dpimap[dpi]
    #print self.dpi

  def checkbox(self,tag,state):
    if tag == "Enable Ray Tracing":
      Nbuttons = self.raymode.numbuttons()
      if state:
        self.ray = 1
        self.RayTrace.button(1).config(state='normal')
        for i in range(Nbuttons):
          self.raymode.button(i).config(state='normal')
      else:
        self.ray = 0
        self.RayTrace.button(1).config(state='disabled')
        for i in range(Nbuttons):
          self.raymode.button(i).config(state='disabled')

    if tag == "Self Shadowing":
      if state:
        cmd.set("ray_shadow", "1")
      else:
        cmd.set("ray_shadow", "0")
    #print self.ray

  def raymode(self,mode):
    if mode == self.modes[0]:
      self.mode = 0
    if mode == self.modes[1]:
      self.mode = 1
    if mode == self.modes[2]:
      self.mode = 3
    cmd.set("ray_trace_mode", str(self.mode))

################################################################################

id_t = "(name,resn,chain,resi)"

def SetUpEnv(obj_name):
  #set ups environment
  cmd.set ("dot_solvent","1")
  cmd.set ("dot_density","4")
  if stored.PDI_Viz_RemoveHETATM == 1:
    cmd.remove(obj_name+" and "+"hetatm")
    cmd.remove(obj_name+" and "+"name MG")
    cmd.remove(obj_name+" and "+"name ZN")
  if stored.PDI_Viz_RemoveH == 1:
    cmd.remove(obj_name+" and "+"h.")
  if stored.PDI_Viz_RemoveAltLocs == 1:
    RemoveAltAtoms(obj_name)
  FixVDW(obj_name)
  FixName(obj_name)

def StoreOriginalAtomID(obj_name):
  stored.PDI_Viz_OrX = {}
  #store original atom ID
  cmd.iterate(obj_name,'stored.PDI_Viz_OrX['+id_t+'] = ID')

def CalculateSplitASA(obj_name):
  stored.PDI_Viz_T_ASA = {}
  stored.PDI_Viz_PROTEIN_ASA = {}
  stored.PDI_Viz_DNA_ASA = {}
  stored.PDI_Viz_DNA_BB_ASA = {}
  stored.PDI_Viz_DNA_BB_PROTEIN_ASA = {}
  stored.PDI_Viz_DNA_BASES_ASA = {}
  stored.PDI_Viz_DNA_BASES_PROTEIN_ASA = {}
  #new mode 3 models
  stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA = {}
  stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA = {}
  stored.PDI_Viz_DNA_MINOR_GROOVE_ASA = {}
  stored.PDI_Viz_DNA_MAJOR_GROOVE_ASA = {}
  #selection names
  dna_bb = stored.PDI_Viz_DNA_bb
  dna_ba = stored.PDI_Viz_DNA_base
  prot_bb = stored.PDI_Viz_p_bb
  prot_sc = stored.PDI_Viz_p_sc
  t_dna = dna_bb+' + '+dna_ba
  t_prot = prot_bb+' + '+prot_sc
  dna_bb_prot = dna_bb+' + '+t_prot
  dna_ba_prot = dna_ba+' + '+t_prot
  #mode 3
  dna_majorg = stored.PDI_Viz_DNA_majorg
  dna_minorg = stored.PDI_Viz_DNA_minorg
  dna_amb = stored.PDI_Viz_DNA_amb

  #dna_no_majorg_prot = dna_minorg+' + '+dna_amb+' + '+t_prot+' + '+dna_bb
  #dna_no_minorg_prot = dna_majorg+' + '+dna_amb+' + '+t_prot+' + '+dna_bb

  dna_no_majorg_prot = dna_minorg+' + '+t_prot+' + '+dna_bb+' + '+dna_amb
  dna_no_minorg_prot = dna_majorg+' + '+t_prot+' + '+dna_bb+' + '+dna_amb

  #       0    1    2     3    4   5
  tpl = "(name,resn,chain,resi,vdw,b)"

  #calculates total area and store original q and b

  #store q and b

  stored.PDI_Viz_or_b = {}
  stored.PDI_Viz_or_q = {}
  cmd.iterate(obj_name,"stored.PDI_Viz_or_b[ID]=b")
  cmd.iterate(obj_name,"stored.PDI_Viz_or_q[ID]=q")

  #calculate area
  cmd.get_area(obj_name,load_b=1)
  cmd.iterate(obj_name,'stored.PDI_Viz_T_ASA['+id_t+']='+tpl)


  #restore q and b
  cmd.alter(obj_name,"b = stored.PDI_Viz_or_b[ID]")
  cmd.alter(obj_name,"q = stored.PDI_Viz_or_q[ID]")

  #creates separated protein and dna objects
  cmd.create(obj_name+"_temp_dna",t_dna)
  cmd.create(obj_name+"_temp_prot",t_prot)

  #calculates area
  cmd.get_area(obj_name+"_temp_dna",load_b=1)
  cmd.get_area(obj_name+"_temp_prot",load_b=1)

  cmd.iterate(obj_name+"_temp_dna",'stored.PDI_Viz_DNA_ASA['+id_t+']='+tpl)
  cmd.iterate(obj_name+"_temp_prot",'stored.PDI_Viz_PROTEIN_ASA['+id_t+']='+tpl)

  cmd.delete(obj_name+"_temp_dna")
  cmd.delete(obj_name+"_temp_prot")

  #creates dna backbone and dna base + prot asa objects and calculates area

  cmd.create('prt_dna_ba_obj',dna_ba_prot)
  cmd.create('dna_bb_obj',dna_bb)

  cmd.get_area('prt_dna_ba_obj',load_b=1)
  cmd.get_area('dna_bb_obj',load_b=1)

  cmd.iterate('prt_dna_ba_obj','stored.PDI_Viz_DNA_BASES_PROTEIN_ASA['+id_t+']='+tpl)
  cmd.iterate('dna_bb_obj','stored.PDI_Viz_DNA_BB_ASA['+id_t+']='+tpl)

  cmd.delete('prt_dna_ba_obj')
  cmd.delete('dna_bb_obj')

  #creates dna bases and dna backbone + prot asa objects and calculates area
  cmd.create('prt_dna_bb_obj',dna_bb_prot)
  cmd.create('dna_ba_obj',dna_ba)

  cmd.get_area('prt_dna_bb_obj',load_b=1)
  cmd.get_area('dna_ba_obj',load_b=1)

  cmd.iterate('prt_dna_bb_obj','stored.PDI_Viz_DNA_BB_PROTEIN_ASA['+id_t+']='+tpl)
  cmd.iterate('dna_ba_obj','stored.PDI_Viz_DNA_BASES_ASA['+id_t+']='+tpl)

  cmd.delete('prt_dna_bb_obj')
  cmd.delete('dna_ba_obj')

  #create dna_groove objects and calculate area

  cmd.create('dna_minorg_',dna_minorg)
  cmd.create('dna_majorg_',dna_majorg)

  cmd.get_area('dna_minorg_',load_b=1)
  cmd.get_area('dna_majorg_',load_b=1)

  cmd.iterate('dna_minorg_','stored.PDI_Viz_DNA_MINOR_GROOVE_ASA['+id_t+']='+tpl)
  cmd.iterate('dna_majorg_','stored.PDI_Viz_DNA_MAJOR_GROOVE_ASA['+id_t+']='+tpl)

  cmd.delete('dna_minorg_')
  cmd.delete('dna_majorg_')

  #create dna_groove + prot objects and calculate area
  cmd.create('dna_no_minorg_prot',dna_no_minorg_prot)
  cmd.create('dna_no_majorg_prot',dna_no_majorg_prot)

  cmd.get_area('dna_no_minorg_prot',load_b=1)
  cmd.get_area('dna_no_majorg_prot',load_b=1)

  cmd.iterate('dna_no_minorg_prot','stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA['+id_t+']='+tpl)
  cmd.iterate('dna_no_majorg_prot','stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA['+id_t+']='+tpl)

  cmd.delete('dna_no_minorg_prot')
  cmd.delete('dna_no_majorg_prot')

  #store results
  stored.PDI_Viz_Results[obj_name] = [stored.PDI_Viz_T_ASA,
                    stored.PDI_Viz_PROTEIN_ASA,
                    stored.PDI_Viz_DNA_ASA,
                    stored.PDI_Viz_DNA_BB_ASA,
                    stored.PDI_Viz_DNA_BB_PROTEIN_ASA,
                    stored.PDI_Viz_DNA_BASES_ASA,
                    stored.PDI_Viz_DNA_BASES_PROTEIN_ASA,
                    stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA,
                    stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA,
                    stored.PDI_Viz_DNA_MAJOR_GROOVE_ASA,
                    stored.PDI_Viz_DNA_MINOR_GROOVE_ASA]


def StoreAreas(obj_name):
  stored.PDI_Viz_Stats[obj_name] = {}
  st = stored.PDI_Viz_Stats[obj_name]

  labels = ['Complex SASA','Free protein SASA','Free DNA SASA',
    'Protein backbone SASA','Protein side chain SASA','DNA backbone SASA','DNA bases SASA',
    'Buried protein surface','Buried protein backbone surface','Buried protein side chain surface',
    'Buried DNA surface','Buried DNA backbone surface','Buried DNA bases surface',
    'DNA major groove surface','DNA minor groove surface',
    'Buried DNA major groove surface','Buried DNA minor groove surface',
    'Buried protein polar surface','Buried protein apolar surface',
    'Buried DNA polar surface','Buried DNA apolar surface',
    'Interface Area']

  keys = ['total','protein','dna',
        'prot_bb','prot_sc', 'dna_bb','dna_ba',
        'prot-total','prot-total_bb','prot-total_sc',
        'dna-total','dna-total_bb','dna-total_ba',
        'dna_ma','dna_mi','dna-total_ma','dna-total_mi',
        'prot-total_prot_polar','prot-total_prot_apolar',
        'dna-total_dna_polar','dna-total_dna_apolar',
        'interface']
  st['total']= "%.1f" % sum(i[5] for i in list(stored.PDI_Viz_Results[obj_name][0].values()))
  st['protein']= "%.1f" % sum(i[5] for i in list(stored.PDI_Viz_Results[obj_name][1].values()))
  st['dna']= "%.1f" % sum(i[5] for i in list(stored.PDI_Viz_Results[obj_name][2].values()))

  #free protein side chains and backbone surfaces
  st['prot_bb']= "%.1f" % sum(stored.PDI_Viz_Results[obj_name][1][item][5] for item in stored.S_IDs[obj_name]['p_bb'])
  st['prot_sc']= "%.1f" % sum(stored.PDI_Viz_Results[obj_name][1][item][5] for item in stored.S_IDs[obj_name]['p_sc'])


  #free dna bases and bb surface
  st['dna_bb']= "%.1f" % sum(stored.PDI_Viz_Results[obj_name][2][item][5] for item in stored.S_IDs[obj_name]['dna_bb'])
  st['dna_ba']= "%.1f" % sum(stored.PDI_Viz_Results[obj_name][2][item][5] for item in stored.S_IDs[obj_name]['dna_ba'])

  #buried protein surface
  bprot = sum(stored.PDI_Viz_PROTEIN_m_TOTAL_d.values())
  st['prot-total'] = "%.1f" % bprot
  #buried protein backbone and base surfaces
  st['prot-total_bb'] = "%.1f" % sum(stored.PDI_Viz_PROTEIN_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['p_bb'])
  st['prot-total_sc'] = "%.1f" % sum(stored.PDI_Viz_PROTEIN_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['p_sc'])
  #polar and apolar areas
  #st['total_prot_apolar'] = "%.1f" % 0.0
  #st['total_prot_polar'] = "%.1f" % 0.0
  st['prot-total_prot_apolar'] = "%.1f" % sum(stored.PDI_Viz_PROTEIN_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['p_apolar'])
  st['prot-total_prot_polar'] = "%.1f" % sum(stored.PDI_Viz_PROTEIN_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['p_polar'])
  
  #buried dna surface
  bdna = sum(stored.PDI_Viz_DNA_m_TOTAL_d.values())
  st['dna-total'] = "%.1f" % bdna
  #buried dna babckbone and base surfaces
  st['dna-total_bb'] = "%.1f" % sum(stored.PDI_Viz_DNA_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['dna_bb'])
  st['dna-total_ba'] = "%.1f" % sum(stored.PDI_Viz_DNA_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['dna_ba'])
  #groove areas
  st['dna_ma']= "%.1f" % sum(stored.PDI_Viz_Results[obj_name][2][item][5] for item in stored.S_IDs[obj_name]['dna_ma'])
  st['dna_mi']= "%.1f" % sum(stored.PDI_Viz_Results[obj_name][2][item][5] for item in stored.S_IDs[obj_name]['dna_mi'])
  #groove dsasa
  st['dna-total_ma'] = "%.1f" % sum(stored.PDI_Viz_DNA_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['dna_ma'])
  st['dna-total_mi'] = "%.1f" % sum(stored.PDI_Viz_DNA_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['dna_mi'])

  #polar and apolar areas
  #st['total_dna_apolar'] = "%.1f" % 0.0
  #st['total_dna_polar'] = "%.1f" % 0.0
  st['dna-total_dna_apolar'] = "%.1f" % sum(stored.PDI_Viz_DNA_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['dna_apolar'])
  st['dna-total_dna_polar'] =  "%.1f" % sum(stored.PDI_Viz_DNA_m_TOTAL_d[item] for item in stored.S_IDs[obj_name]['dna_polar'])

  #interface area
  st['interface'] = "%.1f" % ((bprot + bdna)/2.0)
  for i in range(len(labels)):
    print(labels[i]+":      ",st[keys[i]])
  print(obj_name,"done.")

def PrintAreas(obj_name):
  pass

def CalculateDelta():
  stored.PDI_Viz_DNA_m_TOTAL_d = {}
  stored.PDI_Viz_PROTEIN_m_TOTAL_d = {}
  stored.PDI_Viz_PROTEIN_m_PROTEIN_BB_d = {}
  stored.PDI_Viz_PROTEIN_m_PROTEIN_BA_d = {}
  stored.PDI_Viz_PROTEIN_NOMAJORG_DNA_m_TOTAL_d = {}
  stored.PDI_Viz_PROTEIN_NOMINORG_DNA_m_TOTAL_d = {}
  stored.PDI_Viz_PROTEIN_BASE_m_TOTAL_d = {}
  stored.PDI_Viz_PROTEIN_BB_m_TOTAL_d = {}
  stored.PDI_Viz_BSA = {}
  #calculates delta between DNA and complete structure, meaning the contact surface area in the DNA
  for key in stored.PDI_Viz_DNA_ASA:
    stored.PDI_Viz_DNA_m_TOTAL_d[key] = stored.PDI_Viz_DNA_ASA[key][5] - stored.PDI_Viz_T_ASA[key][5]
  #calculates delta between protein and complete structure, meaning the contact surface area in the Protein
  for key in stored.PDI_Viz_PROTEIN_ASA:
    stored.PDI_Viz_PROTEIN_m_TOTAL_d[key] = stored.PDI_Viz_PROTEIN_ASA[key][5] - stored.PDI_Viz_T_ASA[key][5]
  #save both deltas
  stored.PDI_Viz_BSA.update(stored.PDI_Viz_DNA_m_TOTAL_d)
  stored.PDI_Viz_BSA.update(stored.PDI_Viz_PROTEIN_m_TOTAL_d)
  #delta between protein and protein+dna bb,
  for key in stored.PDI_Viz_PROTEIN_ASA:
    stored.PDI_Viz_PROTEIN_m_PROTEIN_BB_d[key] = stored.PDI_Viz_PROTEIN_ASA[key][5] - stored.PDI_Viz_DNA_BB_PROTEIN_ASA[key][5]
  #delta between protein and protein+dna ba
  for key in stored.PDI_Viz_PROTEIN_ASA:
    stored.PDI_Viz_PROTEIN_m_PROTEIN_BA_d[key] = stored.PDI_Viz_PROTEIN_ASA[key][5] - stored.PDI_Viz_DNA_BASES_PROTEIN_ASA[key][5]
  #delta between protein+dna wo major groove and complex
  for key in stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA:
    stored.PDI_Viz_PROTEIN_NOMAJORG_DNA_m_TOTAL_d[key] = stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA[key][5] - stored.PDI_Viz_T_ASA[key][5]
  #delta between protein+dna wo minor groove and complex
  for key in stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA:
    stored.PDI_Viz_PROTEIN_NOMINORG_DNA_m_TOTAL_d[key] = stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA[key][5] - stored.PDI_Viz_T_ASA[key][5]
  #delta between protein+bases and complex
  for key in stored.PDI_Viz_DNA_BASES_PROTEIN_ASA:
    stored.PDI_Viz_PROTEIN_BASE_m_TOTAL_d[key] = stored.PDI_Viz_DNA_BASES_PROTEIN_ASA[key][5] - stored.PDI_Viz_T_ASA[key][5]
  #delta between protein+backbone and complex
  for key in stored.PDI_Viz_DNA_BB_PROTEIN_ASA:
    stored.PDI_Viz_PROTEIN_BB_m_TOTAL_d[key] = stored.PDI_Viz_DNA_BB_PROTEIN_ASA[key][5] - stored.PDI_Viz_T_ASA[key][5]

def Mode(mode,obj_name,bg,dna_c,dna_bb_c,dna_ba_c,dna_both,T,prot_grad,dna_noint,dna_acceptor,dna_donor,dna_thymine,dna_neutral):
  available_modes = list(range(12))
  #print "Mode ",mode," selected."
  mode = int(mode)
  if(mode in available_modes):
    #create models
    if mode in [0,1,2,3]:
      stored.PDI_Viz_prot_obj = obj_name + "_prot_obj_m2"
      cmd.create(stored.PDI_Viz_prot_obj,stored.PDI_Viz_Prot)
      cmd.alter("object "+stored.PDI_Viz_prot_obj,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj, 'b = stored.PDI_Viz_PROTEIN_m_PROTEIN_BB_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj, 'q = stored.PDI_Viz_PROTEIN_m_PROTEIN_BA_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj,'partial_charge = stored.PDI_Viz_PROTEIN_m_TOTAL_d['+id_t+']')
      cmd.select(stored.PDI_Viz_qGTZeroProt, 'q > '+str(stored.PDI_Viz_dSASAcutoff)+' and object '+stored.PDI_Viz_prot_obj)
      cmd.select(stored.PDI_Viz_bGTZeroProt, 'b > '+str(stored.PDI_Viz_dSASAcutoff)+' and object '+ stored.PDI_Viz_prot_obj)
    if mode in [0,1,2]:
      stored.PDI_Viz_dna_obj = obj_name + "_dna_obj_m3"
      cmd.create(stored.PDI_Viz_dna_obj,stored.PDI_Viz_DNA)
      cmd.alter("object "+stored.PDI_Viz_dna_obj,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.select('dna_bb_obj_','object '+stored.PDI_Viz_dna_obj+' and name '+','.join(stored.PDI_Viz_dna_backbone))
      cmd.select('dna_ba_obj_','object '+stored.PDI_Viz_dna_obj+' and name '+','.join(stored.PDI_Viz_dna_base))
      cmd.alter("object "+stored.PDI_Viz_dna_obj,'b = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj,'q = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj,'partial_charge = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
      cmd.select(stored.PDI_Viz_bGTZeroDNA, 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_bb_obj_')
      cmd.select(stored.PDI_Viz_qGTZeroDNA, 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_ba_obj_')
    if mode in [0,3]:
      stored.PDI_Viz_dna_obj4 = obj_name + "_dna_obj_m2"
      cmd.create(stored.PDI_Viz_dna_obj4,stored.PDI_Viz_DNA)
      cmd.alter("object "+stored.PDI_Viz_dna_obj4,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.select('dna_bb_obj2_','object '+stored.PDI_Viz_dna_obj4+' and name '+','.join(stored.PDI_Viz_dna_backbone))
      cmd.select('dna_ba_obj2_','object '+stored.PDI_Viz_dna_obj4+' and name '+','.join(stored.PDI_Viz_dna_base))
      cmd.alter("object "+stored.PDI_Viz_dna_obj4,'b = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj4,'q = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj4,'partial_charge = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
      cmd.select("bGTZeroDNA_4", 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_bb_obj2_')
      cmd.select("qGTZeroDNA_4", 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_ba_obj2_')

    if mode in [0,4,5,6]:
      stored.PDI_Viz_prot_obj3 = obj_name + "_prot_obj_g3" #protein object for major/minor groove
      cmd.create(stored.PDI_Viz_prot_obj3,stored.PDI_Viz_Prot)
      cmd.alter("object "+stored.PDI_Viz_prot_obj3,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj3, 'q = stored.PDI_Viz_PROTEIN_NOMAJORG_DNA_m_TOTAL_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj3, 'b = stored.PDI_Viz_PROTEIN_NOMINORG_DNA_m_TOTAL_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj3, 'partial_charge = stored.PDI_Viz_PROTEIN_m_TOTAL_d['+id_t+']')
      cmd.select(stored.PDI_Viz_bGTZeroProt2, 'b > '+str(stored.PDI_Viz_dSASAcutoff)+' and object '+ stored.PDI_Viz_prot_obj3)
      cmd.select(stored.PDI_Viz_qGTZeroProt2, 'q > '+str(stored.PDI_Viz_dSASAcutoff)+' and object '+stored.PDI_Viz_prot_obj3)

    if mode in [0,4,5]:
      stored.PDI_Viz_dna_obj3 = obj_name + "_dna_obj_g3"   #dna object for major/minor groove
      cmd.create(stored.PDI_Viz_dna_obj3,stored.PDI_Viz_DNA)
      cmd.alter("object "+stored.PDI_Viz_dna_obj3,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.select('dna_majorg',GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj3,stored.PDI_Viz_dna_majorgroove))
      cmd.select('dna_minorg',GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj3,stored.PDI_Viz_dna_minorgroove))
      cmd.alter("object "+stored.PDI_Viz_dna_obj3,'b = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj3,'q = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj3,'partial_charge = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
      cmd.select(stored.PDI_Viz_bGTZeroDNA2, 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_minorg')
      cmd.select(stored.PDI_Viz_qGTZeroDNA2, 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_majorg')


    if mode in [0,7,8,9]:
      stored.PDI_Viz_prot_obj2 = obj_name + "_prot_int_m2" #chemical interaction
      cmd.create(stored.PDI_Viz_prot_obj2,stored.PDI_Viz_Prot)
      cmd.alter("object "+stored.PDI_Viz_prot_obj2,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj2, 'b = stored.PDI_Viz_PROTEIN_m_PROTEIN_BB_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj2, 'q = stored.PDI_Viz_PROTEIN_m_PROTEIN_BA_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj2, 'partial_charge = stored.PDI_Viz_PROTEIN_m_TOTAL_d['+id_t+']')
      prot_HA = '('+GenSeleFromResnNameDict(stored.PDI_Viz_prot_obj2,stored.PDI_Viz_prot_Hbond_A)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      prot_HD = '('+GenSeleFromResnNameDict(stored.PDI_Viz_prot_obj2,stored.PDI_Viz_prot_Hbond_D)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      prot_Hn = '(object '+stored.PDI_Viz_prot_obj2+ ') and (not ('+prot_HA+')) and (not ('+prot_HD+')) and not (name '+','.join(stored.PDI_Viz_prot_backbone)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      prot_chemsel = ["prot_HA","prot_HD","prot_Hn"]
      cmd.select(prot_chemsel[0],   prot_HA)
      cmd.select(prot_chemsel[1],   prot_HD)
      cmd.select(prot_chemsel[2],   prot_Hn)
    if mode in [0,7,8]:
      stored.PDI_Viz_dna_obj2 =obj_name + "_dna_int_m2"    #chemical interaction
      cmd.create(stored.PDI_Viz_dna_obj2,stored.PDI_Viz_DNA)
      cmd.alter("object "+stored.PDI_Viz_dna_obj2,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_dna_obj2,'b = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj2,'q = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj2,'partial_charge = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
      dna_HA = '('+GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj2,stored.PDI_Viz_dna_Hbond_A)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      dna_HD = '('+GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj2,stored.PDI_Viz_dna_Hbond_D)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      dna_Tmet = stored.PDI_Viz_dna_obj2 + ' and resn DT and name C7 and partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)
      dna_Hn = '('+GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj2,stored.PDI_Viz_dna_Hn)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      dna_chemsel = ["dna_HA","dna_HD","dna_Tmet", "dna_Hn"]
      cmd.select(dna_chemsel[0],   dna_HA)
      cmd.select(dna_chemsel[1],   dna_HD)
      cmd.select(dna_chemsel[2], dna_Tmet)
      cmd.select(dna_chemsel[3],   dna_Hn)
    #residue colored objects
    if mode in [0,6]:
      stored.PDI_Viz_dna_obj6 = obj_name + "_dna_obj_g4"
      cmd.create(stored.PDI_Viz_dna_obj6,stored.PDI_Viz_DNA)
      cmd.alter("object "+stored.PDI_Viz_dna_obj6,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.select('dna_majorg2',GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj6,stored.PDI_Viz_dna_majorgroove))
      cmd.select('dna_minorg2',GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj6,stored.PDI_Viz_dna_minorgroove))
      cmd.alter("object "+stored.PDI_Viz_dna_obj6,'b = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj6,'q = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj6,'partial_charge = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
      cmd.select(stored.PDI_Viz_bGTZeroDNA3, 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_minorg2')
      cmd.select(stored.PDI_Viz_qGTZeroDNA3, 'partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+' and dna_majorg2')
    if mode in [0,9]:
      stored.PDI_Viz_dna_obj5 = obj_name + "_dna_int_m3"
      cmd.create(stored.PDI_Viz_dna_obj5,stored.PDI_Viz_DNA)
      cmd.alter("object "+stored.PDI_Viz_dna_obj5,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_dna_obj5,'b = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj5,'q = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj5,'partial_charge = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
      dna_HA_s = '('+GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj5,stored.PDI_Viz_dna_Hbond_A)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      dna_HD_s = '('+GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj5,stored.PDI_Viz_dna_Hbond_D)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      dna_Tmet_s = stored.PDI_Viz_dna_obj5 + ' and resn DT and name C7 and partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)
      dna_Hn_s = '('+GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj5,stored.PDI_Viz_dna_Hn)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      dna_chemsel_seq = ["dna_HA_s","dna_HD_s","dna_Tmet_s", "dna_Hn_s"]
      cmd.select(dna_chemsel_seq[0],   dna_HA_s)
      cmd.select(dna_chemsel_seq[1],   dna_HD_s)
      cmd.select(dna_chemsel_seq[2], dna_Tmet_s)
      cmd.select(dna_chemsel_seq[3],   dna_Hn_s)

    if mode in [11]:
      stored.PDI_Viz_dna_obj_polar = obj_name+ "_dna_polar"
      cmd.create(stored.PDI_Viz_dna_obj_polar,stored.PDI_Viz_DNA)
      cmd.alter("object "+stored.PDI_Viz_dna_obj_polar,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_dna_obj_polar,'b = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj_polar,'q = 0.0')
      cmd.alter("object "+stored.PDI_Viz_dna_obj_polar,'partial_charge = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
      dna_polar_ba = GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj_polar,stored.PDI_Viz_dna_polar_ba)
      dna_polar_bb = GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj_polar,stored.PDI_Viz_dna_polar_bb)
      dna_polar = '(('+dna_polar_ba+') or ('+dna_polar_bb+')) and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      dna_apolar_ba = GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj_polar,stored.PDI_Viz_dna_apolar_ba)
      dna_apolar_bb = GenSeleFromResnNameDict(stored.PDI_Viz_dna_obj_polar,stored.PDI_Viz_dna_apolar_bb)
      dna_apolar = '(('+dna_apolar_ba+') or ('+dna_apolar_bb+')) and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')'
      cmd.select("dna_polar_",dna_polar)
      cmd.select("dna_apolar_",dna_apolar)
      
      stored.PDI_Viz_prot_obj_polar = obj_name + "_prot_polar"
      cmd.create(stored.PDI_Viz_prot_obj_polar,stored.PDI_Viz_Prot)
      cmd.alter("object "+stored.PDI_Viz_prot_obj_polar,'ID = stored.PDI_Viz_OrX['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj_polar, 'b = stored.PDI_Viz_PROTEIN_m_PROTEIN_BB_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj_polar, 'q = stored.PDI_Viz_PROTEIN_m_PROTEIN_BA_d['+id_t+']')
      cmd.alter("object "+stored.PDI_Viz_prot_obj_polar,'partial_charge = stored.PDI_Viz_PROTEIN_m_TOTAL_d['+id_t+']')
      cmd.select("prot_polar_",  '('+GenSeleFromResnNameDict(stored.PDI_Viz_prot_obj_polar,stored.PDI_Viz_prot_polar)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')')
      cmd.select("prot_apolar_", '('+GenSeleFromResnNameDict(stored.PDI_Viz_prot_obj_polar,stored.PDI_Viz_prot_apolar)+') and (partial_charge > '+str(stored.PDI_Viz_dSASAcutoff)+')')
      
    cmd.disable(stored.PDI_Viz_qGTZeroDNA)

    cmd.delete('dna_bb_obj_')
    cmd.delete('dna_ba_obj_')

    cmd.delete('dna_bb_obj2_')
    cmd.delete('dna_ba_obj2_')

    cmd.delete('dna_minorg')
    cmd.delete('dna_majorg')

    cmd.delete('dna_minorg2')
    cmd.delete('dna_majorg2')
    ######## dna color block

    Colorize(obj_name,mode,bg,dna_c,dna_bb_c,dna_ba_c,dna_both,T,prot_grad,dna_noint,dna_acceptor,dna_donor,dna_thymine,dna_neutral)
  else:
    print("Invalid mode.",mode)

def ColorOut(obj_name,path='./',prot = 0,d_byres = 1,d_atom = 0):
  #create models
  PDI_Viz_prot_obj = obj_name + "_prot_obj_t"
  PDI_Viz_dna_obj = obj_name + "_dna_obj_t"

  qGTZeroDNA = obj_name+'_qGTZeroDNA'
  bGTZeroDNA = obj_name+'_bGTZeroDNA'
  qGTZeroProt = obj_name+'_qGTZeroProt'
  bGTZeroProt = obj_name+'_bGTZeroProt'

  cmd.create(PDI_Viz_dna_obj,stored.PDI_Viz_DNA)
  cmd.create(PDI_Viz_prot_obj,stored.PDI_Viz_Prot)
  #select atoms in objects
  cmd.select('dna_bb_obj_t',PDI_Viz_dna_obj + ' and name '+','.join(stored.PDI_Viz_dna_backbone))
  cmd.select('dna_ba_obj_t',PDI_Viz_dna_obj + ' and name '+','.join(stored.PDI_Viz_dna_base))
  #modify b and q in the new objects
  cmd.alter('dna_bb_obj_t','b = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')
  cmd.alter('dna_bb_obj_t','q = 0.0')

  cmd.alter('dna_ba_obj_t','b = 0.0')
  cmd.alter('dna_ba_obj_t','q = stored.PDI_Viz_DNA_m_TOTAL_d['+id_t+']')

  cmd.delete('dna_bb_obj_t')
  cmd.delete('dna_ba_obj_t')

  cmd.alter(PDI_Viz_prot_obj, 'b = stored.PDI_Viz_PROTEIN_m_PROTEIN_BB_d['+id_t+']')
  cmd.alter(PDI_Viz_prot_obj, 'q = stored.PDI_Viz_PROTEIN_m_PROTEIN_BA_d['+id_t+']')

  cmd.alter(PDI_Viz_prot_obj,'partial_charge = b+q')
  cmd.alter(PDI_Viz_dna_obj,'partial_charge = b+q')


  cmd.select(qGTZeroProt, 'q > '+str(stored.PDI_Viz_dSASAcutoff)+' and '+ "object "+PDI_Viz_prot_obj)
  cmd.select(bGTZeroProt, 'b > '+str(stored.PDI_Viz_dSASAcutoff)+' and '+ "object "+PDI_Viz_prot_obj)

  cmd.select(bGTZeroDNA, 'b > '+str(stored.PDI_Viz_dSASAcutoff)+' and '+ PDI_Viz_dna_obj)
  cmd.select(qGTZeroDNA, 'q > '+str(stored.PDI_Viz_dSASAcutoff)+' and '+ PDI_Viz_dna_obj)


  ns = {'dna':{},'prot':{}}
  cmd.iterate('byres '+stored.PDI_Viz_DNA,"dna[(name,resn,chain,resi)]='NA'",space=ns)
  cmd.iterate('byres '+bGTZeroDNA,"dna[(name,resn,chain,resi)]='BB'",space=ns)
  cmd.iterate('byres '+qGTZeroDNA,"dna[(name,resn,chain,resi)]='BASE'",space=ns)
  cmd.iterate('byres '+bGTZeroDNA+' and '+'byres '+qGTZeroDNA,"dna[(name,resn,chain,resi)]='BOTH'",space=ns)

  cmd.iterate(stored.PDI_Viz_Prot,"prot[(name,resn,chain,resi)]='NA'",space=ns)
  cmd.iterate(bGTZeroProt,"prot[(name,resn,chain,resi)]='BB'",space=ns)
  cmd.iterate(qGTZeroProt,"prot[(name,resn,chain,resi)]='BASE'",space=ns)
  cmd.iterate(bGTZeroProt+' and '+qGTZeroProt,"prot[(name,resn,chain,resi)]='BOTH'",space=ns)

  if (d_byres != 0):
    dnafile = open(path+obj_name+'_pdiviz_dna_babb_color_byres.csv','w')
    dna = []
    for key in ns['dna']:
      dna.append([key[0],key[1],key[2],key[3],ns['dna'][key]])
    dna.sort(key=lambda a: int(a[3]))
    temp_dna_list = []
    d_res_list = []
    for item in dna:
      if(not [item[2],item[3]] in d_res_list):
        d_res_list.append([item[2],item[3]])
        temp_dna_list.append([item[1],item[2],item[3],item[4]])
    for item in temp_dna_list:
      for s in item:
        dnafile.write(str(s)+"\t")
      dnafile.write(Lend)
    dnafile.close()

  cmd.iterate(stored.PDI_Viz_DNA,"dna[(name,resn,chain,resi)]='n'",space=ns)
  cmd.iterate(bGTZeroDNA,"dna[(name,resn,chain,resi)]='bb'",space=ns)
  cmd.iterate(qGTZeroDNA,"dna[(name,resn,chain,resi)]='ba'",space=ns)
  cmd.iterate(bGTZeroDNA+' and '+qGTZeroDNA,"dna[(name,resn,chain,resi)]='both'",space=ns)

  if(d_atom != 0):
    dnafile = open(path+obj_name+'_pdiviz_dna_babb_color_byatom.csv','w')
    dna = []
    for key in ns['dna']:
      dna.append([key[0],key[1],key[2],key[3],ns['dna'][key]])
    dna.sort(key=lambda a: a[3])
    for item in dna:
      for s in item:
        dnafile.write(str(s)+"\t")
      dnafile.write(Lend)
    dnafile.close()

  if(prot != 0):
    prot=[]
    protfile = open(path+obj_name+'_pdiviz_prot_babb_color_byatom.csv','w')
    for key in ns['prot']:
      prot.append([key[0],key[1],key[2],key[3],ns['prot'][key]])
    prot.sort(key=lambda a: a[3])
    for item in prot:
      for s in item:
        protfile.write(str(s)+"\t")
      protfile.write(Lend)
    protfile.close()
  #destroy models
  cmd.delete(PDI_Viz_prot_obj)
  cmd.delete(PDI_Viz_dna_obj)
  cmd.delete(bGTZeroDNA)
  cmd.delete(qGTZeroDNA)
  cmd.delete(bGTZeroProt)
  cmd.delete(qGTZeroProt)


def Colorize(obj_name,mode,bg,dna_c,dna_bb_c,dna_ba_c,dna_both,T,prot_grad,dna_noint,dna_acceptor,dna_donor,dna_thymine,dna_neutral):
  mode = int(mode)
  cmd.set('cartoon_ladder_mode', 0)
  cmd.set('sphere_scale', 0.3)
  cmd.bg_color(bg)
  cmd.set('transparency', T)
  #color list and selection for chemical bond mode
  dna_scolorlist = [dna_acceptor,dna_donor,dna_thymine,dna_neutral]
  dna_chemsel = ["dna_HA","dna_HD","dna_Tmet", "dna_Hn"]
  dna_chemsel_s = ["dna_HA_s","dna_HD_s","dna_Tmet_s", "dna_Hn_s"]
  prot_scolorlist = [dna_acceptor,dna_donor,dna_neutral]
  prot_chemsel = ["prot_HA","prot_HD", "prot_Hn"]

  ##DNA
  if mode in [0,1,2]:
    cmd.hide('everything', 'object '+stored.PDI_Viz_dna_obj)
    cmd.show('sticks', 'byres ('+stored.PDI_Viz_bGTZeroDNA+' or '+stored.PDI_Viz_qGTZeroDNA+')')
    cmd.set_bond('stick_radius', 0.1, 'object '+stored.PDI_Viz_dna_obj)
    cmd.show('spheres', '('+stored.PDI_Viz_bGTZeroDNA+' or '+stored.PDI_Viz_qGTZeroDNA+')')
    if mode == 2:
      cmd.show('surface', 'object '+stored.PDI_Viz_dna_obj)
    else:
      cmd.show('cartoon', 'object '+stored.PDI_Viz_dna_obj)
    cmd.color(dna_c, 'object '+stored.PDI_Viz_dna_obj)
    color_b(selection=stored.PDI_Viz_bGTZeroDNA, item='partial_charge', mode='ramp', gradient='user',user_rgb=prot_grad[1], nbins=10, sat=0.7, value=0.7,debug=0)
    color_b(selection=stored.PDI_Viz_qGTZeroDNA, item='partial_charge', mode='ramp', gradient='user',user_rgb=prot_grad[0], nbins=10, sat=0.7, value=0.8,debug=0)
    color_b(selection=(stored.PDI_Viz_bGTZeroDNA + ' and ' + stored.PDI_Viz_qGTZeroDNA), item='partial_charge', mode='ramp', gradient='user', user_rgb=prot_grad[2], nbins=10, sat=0.7, value=0.9,debug=0)
  if mode in [0,7,8]:
    cmd.hide('everything', 'object '+stored.PDI_Viz_dna_obj2)

    cmd.show('sticks', 'byres (object '+stored.PDI_Viz_dna_obj2+' and partial_charge > 0.0)')
    cmd.set_bond('stick_radius', 0.1, 'object '+stored.PDI_Viz_dna_obj2)
    cmd.show('spheres',stored.PDI_Viz_dna_obj2+" and (partial_charge > "+str(stored.PDI_Viz_dSASAcutoff)+") and ("+" or ".join(dna_chemsel)+")")
    if mode == 8:
      cmd.show('surface', 'object '+stored.PDI_Viz_dna_obj2)
    else:
      cmd.show('cartoon', 'object '+stored.PDI_Viz_dna_obj2)
    cmd.color(dna_noint,'object '+stored.PDI_Viz_dna_obj2)
    for item in zip(dna_chemsel,dna_scolorlist):
      cmd.color(item[1], item[0])
    for item in dna_chemsel:
      cmd.delete(item)
  if mode in [0,4,5]:
    cmd.hide('everything','object '+stored.PDI_Viz_dna_obj3)

    cmd.show('sticks', 'byres ('+stored.PDI_Viz_bGTZeroDNA2+' or '+stored.PDI_Viz_qGTZeroDNA2+')')
    cmd.set_bond('stick_radius', 0.1, 'object '+stored.PDI_Viz_dna_obj3)
    cmd.show('spheres', '('+stored.PDI_Viz_bGTZeroDNA2+' or '+stored.PDI_Viz_qGTZeroDNA2+')')
    if mode == 5:
      cmd.show('surface', 'object '+stored.PDI_Viz_dna_obj3)
    else:
      cmd.show('cartoon', 'object '+stored.PDI_Viz_dna_obj3)
    cmd.color(dna_c,'object '+stored.PDI_Viz_dna_obj3)
    color_b(selection=stored.PDI_Viz_bGTZeroDNA2, item='partial_charge', mode='ramp', gradient='wg', nbins=10, sat=0.7, value=0.7,debug=0)
    color_b(selection=stored.PDI_Viz_qGTZeroDNA2, item='partial_charge', mode='ramp', gradient='wb', nbins=10, sat=0.7, value=0.8,debug=0)
    color_b(selection=(stored.PDI_Viz_bGTZeroDNA2 + ' and ' + stored.PDI_Viz_qGTZeroDNA2), item='partial_charge', mode='ramp', gradient='user', user_rgb=prot_grad[2], nbins=10, sat=0.7, value=0.9,debug=0)
  if mode in [0,3]:
    cmd.hide('everything', stored.PDI_Viz_dna_obj4)
    cmd.show('cartoon', stored.PDI_Viz_dna_obj4)
    cmd.show('sticks', 'byres (bGTZeroDNA_4 or qGTZeroDNA_4)')
    cmd.set_bond('stick_radius', 0.1, 'object '+stored.PDI_Viz_dna_obj4)
    cmd.show('spheres', '(bGTZeroDNA_4 or qGTZeroDNA_4)')
    cmd.color(dna_c,'object '+stored.PDI_Viz_dna_obj4)
    cmd.color(dna_bb_c, 'byres bGTZeroDNA_4')
    cmd.color(dna_ba_c, 'byres qGTZeroDNA_4')
    cmd.color(dna_both, '(byres bGTZeroDNA_4) and (byres qGTZeroDNA_4)')
    cmd.delete('bGTZeroDNA_4')
    cmd.delete('qGTZeroDNA_4')
  if mode in [0,9]:
    cmd.hide('everything', 'object '+stored.PDI_Viz_dna_obj5)
    cmd.show('cartoon', 'object '+stored.PDI_Viz_dna_obj5)
    cmd.show('sticks', 'byres (object '+stored.PDI_Viz_dna_obj5+' and partial_charge > 0.0)')
    cmd.set_bond('stick_radius', 0.1, 'object '+stored.PDI_Viz_dna_obj5)
    cmd.show('spheres',stored.PDI_Viz_dna_obj5+" and (partial_charge > "+str(stored.PDI_Viz_dSASAcutoff)+") and ("+" or ".join(dna_chemsel_s)+")")
    cmd.color(dna_noint,'object '+stored.PDI_Viz_dna_obj5)
    cmd.color(dna_acceptor,'byres DNA_HA_s')
    cmd.color(dna_donor,   'byres DNA_HD_s')
    cmd.color(dna_thymine, 'byres DNA_Tmet_s')
    cmd.color(dna_neutral, 'byres DNA_HN_s')
    cmd.color('dna_DA',  '(byres DNA_HD_s) and (byres DNA_HA_s)')
    cmd.color('dna_DH',  '(byres DNA_HD_s) and (byres DNA_Hn_s)')
    cmd.color('dna_AH',  '(byres DNA_HA_s) and (byres DNA_Hn_s)')
    cmd.color('dna_DAH', '(byres DNA_HA_s) and (byres DNA_HD_s) and (byres DNA_Hn_s)')
    cmd.color('dna_DT',  '(byres DNA_HD_s) and (byres DNA_Tmet_s)')
    cmd.color('dna_AT',  '(byres DNA_HA_s) and (byres DNA_Tmet_s)')
    cmd.color('dna_HT',  '(byres DNA_Hn_s) and (byres DNA_Tmet_s)')
    cmd.color('dna_DAT', '(byres DNA_HD_s) and (byres DNA_HA_s) and (byres DNA_Tmet_s)')
    cmd.color('dna_DTH', '(byres DNA_HD_s) and (byres DNA_Tmet_s) and (byres DNA_Hn_s)')
    cmd.color('dna_ATH', '(byres DNA_HA_s) and (byres DNA_Tmet_s) and (byres DNA_Hn_s)')
    cmd.color('dna_DATH','(byres DNA_HD_s) and (byres DNA_HA_s) and (byres DNA_Tmet_s) and (byres DNA_Hn_s)')
    for item in dna_chemsel_s:
      cmd.delete(item)

  if mode in [0,6]:
    cmd.hide('everything','object '+stored.PDI_Viz_dna_obj6)
    cmd.show('cartoon', 'object '+stored.PDI_Viz_dna_obj6)
    cmd.show('sticks', 'byres ('+stored.PDI_Viz_bGTZeroDNA3+' or '+stored.PDI_Viz_qGTZeroDNA3+')')
    cmd.set_bond('stick_radius', 0.1, 'object '+stored.PDI_Viz_dna_obj6)
    cmd.show('spheres', '('+stored.PDI_Viz_bGTZeroDNA3+' or '+stored.PDI_Viz_qGTZeroDNA3+')')
    cmd.color(dna_c,'object '+stored.PDI_Viz_dna_obj6)
    cmd.color('green', 'byres '+stored.PDI_Viz_bGTZeroDNA3)
    cmd.color(dna_ba_c, 'byres '+stored.PDI_Viz_qGTZeroDNA3,)
    cmd.color('purple', '(byres '+stored.PDI_Viz_bGTZeroDNA3+') and (byres '+stored.PDI_Viz_qGTZeroDNA3+')')
    cmd.delete('bGZeroDNA_g4')
    cmd.delete('qGZeroDNA_g4')

  #proteins
  #show sticks, byres bGreaterZeroProtein and not name N+C+O for obj1 (base backbone)
  if mode in [0,1,2,3]:
    cmd.color('white', "object "+stored.PDI_Viz_prot_obj)
    cmd.hide('lines', "object "+stored.PDI_Viz_prot_obj)
    cmd.show('sticks', 'byres ('+stored.PDI_Viz_bGTZeroProt+' or '+stored.PDI_Viz_qGTZeroProt+')')
    cmd.show('spheres', stored.PDI_Viz_bGTZeroProt+' or '+stored.PDI_Viz_qGTZeroProt)
    cmd.set_bond('stick_radius', 0.1, stored.PDI_Viz_prot_obj)
    if mode !=2 :
      cmd.show('surface', "object "+stored.PDI_Viz_prot_obj)
    cmd.show('cartoon', "object "+stored.PDI_Viz_prot_obj)
    color_b(selection=stored.PDI_Viz_bGTZeroProt, item='b', mode='ramp', gradient='user',user_rgb=prot_grad[1], nbins=10, sat=0.7, value=0.7,debug=0)
    color_b(selection=stored.PDI_Viz_qGTZeroProt, item='q', mode='ramp', gradient='user',user_rgb=prot_grad[0], nbins=10, sat=0.7, value=0.8,debug=0)
    color_b(selection=(stored.PDI_Viz_bGTZeroProt + ' and ' + stored.PDI_Viz_qGTZeroProt), item='partial_charge', mode='ramp', gradient='user', user_rgb=prot_grad[2], nbins=10, sat=0.7, value=0.9,debug=0)
  #show sticks, byres bGreaterZeroProtein and not name N+C+O for obj2 (interaction)
  if mode in [0,7,8,9]:
    cmd.color(dna_noint, "object "+stored.PDI_Viz_prot_obj2)
    cmd.hide('lines', "object "+stored.PDI_Viz_prot_obj2)
    cmd.show('sticks', 'byres (object '+stored.PDI_Viz_prot_obj2+ ' and partial_charge > 0)')
    cmd.show('spheres',stored.PDI_Viz_prot_obj2+" and (partial_charge > "+str(stored.PDI_Viz_dSASAcutoff)+") and ("+" or ".join(prot_chemsel)+")")
    cmd.set_bond('stick_radius', 0.1, stored.PDI_Viz_prot_obj2)
    if mode != 8:
      cmd.show('surface', "object "+stored.PDI_Viz_prot_obj2)
    cmd.show('cartoon', "object "+stored.PDI_Viz_prot_obj2)
    for item in zip(prot_chemsel,prot_scolorlist):
      cmd.color(item[1], item[0])
    cmd.color("violet", "prot_HA and prot_HD")
    for item in prot_chemsel:
      cmd.delete(item)

  #show sticks, byres bGreaterZeroProtein and not name N+C+O for obj3 (major minor groove)
  if mode in [0,4,5,6]:
    cmd.color('white', "object "+stored.PDI_Viz_prot_obj3)
    cmd.hide('lines', "object "+stored.PDI_Viz_prot_obj3)
    cmd.show('sticks', 'byres ('+stored.PDI_Viz_bGTZeroProt2+' or '+stored.PDI_Viz_qGTZeroProt2+')')
    cmd.show('spheres', stored.PDI_Viz_bGTZeroProt2+' or '+stored.PDI_Viz_qGTZeroProt2+" and (partial_charge > "+str(stored.PDI_Viz_dSASAcutoff)+")")
    cmd.set_bond('stick_radius', 0.1, stored.PDI_Viz_prot_obj3)
    if mode != 5:
      cmd.show('surface', "object "+stored.PDI_Viz_prot_obj3)
    cmd.show('cartoon', "object "+stored.PDI_Viz_prot_obj3)
    color_b(selection=stored.PDI_Viz_bGTZeroProt2, item='partial_charge', mode='ramp', gradient='wg', nbins=10, sat=0.7, value=0.7,debug=0)
    color_b(selection=stored.PDI_Viz_qGTZeroProt2, item='partial_charge', mode='ramp', gradient='wb', nbins=10, sat=0.7, value=0.8,debug=0)
    color_b(selection=(stored.PDI_Viz_bGTZeroProt2 + ' and ' + stored.PDI_Viz_qGTZeroProt2), item='partial_charge', mode='ramp', gradient='user', user_rgb='(1.0,1.0,1.0, 0.875,0.5,0.875, 0.75,0.00,0.75)', nbins=10, sat=0.7, value=0.9,debug=0)

  if mode in [11]:
    #protein
    cmd.color('white', "object "+stored.PDI_Viz_prot_obj_polar)
    cmd.hide('lines', "object "+stored.PDI_Viz_prot_obj_polar)
    cmd.show('sticks', 'byres ( prot_polar_ or prot_apolar_ )')
    cmd.show('spheres',"(prot_polar_ or prot_apolar_) and (partial_charge > "+str(stored.PDI_Viz_dSASAcutoff)+")")
    cmd.set_bond('stick_radius', 0.1, stored.PDI_Viz_prot_obj_polar)
    cmd.show('surface', "object "+stored.PDI_Viz_prot_obj_polar)
    cmd.show('cartoon', "object "+stored.PDI_Viz_prot_obj_polar)
    color_b(selection="prot_polar_", item='partial_charge', mode='ramp', gradient='wr', nbins=10, sat=0.7, value=0.7,debug=0)
    color_b(selection="prot_apolar_", item='partial_charge', mode='ramp', gradient='wb', nbins=10, sat=0.7, value=0.8,debug=0)
    #dna
    cmd.hide('everything','object '+stored.PDI_Viz_dna_obj_polar)

    cmd.show('sticks', 'byres (dna_polar_ or dna_apolar_)')
    cmd.set_bond('stick_radius', 0.1, 'object '+stored.PDI_Viz_dna_obj_polar)
    cmd.show('spheres', '(dna_polar_ or dna_apolar_)')
    cmd.show('surface', 'object '+stored.PDI_Viz_dna_obj_polar)
    cmd.show('cartoon', 'object '+stored.PDI_Viz_dna_obj_polar)
    cmd.color(dna_c,'object '+stored.PDI_Viz_dna_obj_polar)
    color_b(selection="dna_polar_", item='partial_charge', mode='ramp', gradient='wr', nbins=10, sat=0.7, value=0.7,debug=0)
    color_b(selection="dna_apolar_", item='partial_charge', mode='ramp', gradient='wb', nbins=10, sat=0.7, value=0.8,debug=0)
    cmd.delete("prot_polar_")
    cmd.delete("prot_apolar_")
    cmd.delete("dna_polar_")
    cmd.delete("dna_apolar_")

def Show(obj_name):
  cmd.disable(obj_name)

def WriteSplitASA(obj_name, path = "./"):

  aa_backbone = stored.PDI_Viz_prot_backbone
  aa_sidechain = stored.PDI_Viz_prot_base
  dna_backbone = stored.PDI_Viz_dna_backbone
  dna_majorgroove = stored.PDI_Viz_dna_majorgroove
  dna_minorgroove = stored.PDI_Viz_dna_minorgroove
  dna_ambiguous = stored.PDI_Viz_dna_ambi

  amino_acids = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

  bases = [ 'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'T','U' ] # DNA, RNA

  deoxybase_lookup = { 'A': 'DA', 'C': 'DC', 'G': 'DG', 'T': 'DT' }

  aa_polar_backbone = ['N', 'O', 'OXT']
  aa_hyd_backbone = ['CA', 'C']

  aa_polar_sidechain = ['ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ', 'NH1', 'NH2',
                        'OD1', 'OD2', 'OG', 'OG1', 'OG2', 'OE1', 'OE2', 'OH']
  aa_hyd_sidechain = ['CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3',
                    'CZ', 'CZ2', 'CZ3', 'CH2']

  dna_polar_backbone = ['O3\'', 'O4\'', 'O5\'', 'OP1', 'OP2', 'OP3', 'O2\''] # O2' of RNA
  dna_hyd_backbone = ['C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'']

  dna_polar_majorgroove = {'DA':['N6', 'N7'],
                         'DC':['N4'],
                         'DG':['N7', 'O6'],
                         'DT':['O4'],
                         'U':['O4']} # RNA
  dna_polar_minorgroove = {'DA':['N1', 'N3'],
                         'DC':['O2'],
                         'DG':['N2', 'N3'],
                         'DT':['N3', 'O2'],
                         'U':['N3', 'O2']}
  dna_polar_ambigous = {'DA':['N9'],
                      'DC':['N1', 'N3'],
                      'DG':['N1', 'N9'],
                      'DT':['N1'],
                      'U':['N1']} # RNA

  dna_hyd_majorgroove = {'DA':['C5', 'C6', 'C8'],
                       'DC':['C4', 'C5', 'C6'],
                       'DG':['C5', 'C6', 'C8'],
                       'DT':['C4', 'C5', 'C6', 'C7','C5M'],
                       'U':['C4', 'C5', 'C6']} # RNA
  dna_hyd_minorgroove = {'DA':['C2', 'C4'],
                       'DC':['C2'],
                       'DG':['C2', 'C4'],
                       'DT':['C2'],
                       'U':['C2']} # RNA
  dna_hyd_ambigous = {'DA':[],
                    'DC':[],
                    'DG':[],
                    'DT':[],
                    'U':[]} # RNA




  formatOutput = lambda resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA : '%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % \
        (resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[0], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[1], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[2], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[3], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[4], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[5], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[6], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[7], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[8],
         resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[9], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[10], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[11], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[12], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[13], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[14],
         resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[15], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[16], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[17], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[18], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[19], resname_chain_resnum_total_ASA_bb_ASA_sc_ASA_majorgroove_ASA_minorgroove_ASA_nogroove_ASA_polar_ASA_polar_bb_ASA_polar_sc_ASA_polar_majorgroove_ASA_polar_minorgroove_ASA_polar_nogroove_ASA_hyd_ASA_hyd_bb_ASA_hyd_sc_ASA_hyd_majorgroove_ASA_hyd_minorgroove_ASA_hyd_nogroove_ASA[20])

  stored.xyz = {}
  cmd.iterate_state(1,"object "+obj_name,"stored.xyz[ID] = (x,y,z)")

  asa_l = [(path+obj_name+"_pdiviz.asa",                                path+obj_name+"_pdiviz.atmasa",                  list(stored.PDI_Viz_T_ASA.values()), 0),
           (path+obj_name+"_pdiviz_protein.asa",                        path+obj_name+"_pdiviz_protein.atmasa",          list(stored.PDI_Viz_PROTEIN_ASA.values()), 0),
           (path+obj_name+"_pdiviz_dna.asa",                            path+obj_name+"_pdiviz_dna.atmasa",              list(stored.PDI_Viz_DNA_ASA.values()), 0),
           (path+obj_name+"_pdiviz_bsa.asa",                             path+obj_name+"_pdiviz_bsa.atmasa",               stored.PDI_Viz_BSA, 1),
#           (path+obj_name+"_pdiviz_dna_bb.asa",                         path+obj_name+"_pdiviz_dna_bb.atmasa",           stored.PDI_Viz_DNA_BB_ASA.values()),
#           (path+obj_name+"_pdiviz_dna_bb_protein.asa",                 path+obj_name+"_pdiviz_dna_bb_protein.atmasa",   stored.PDI_Viz_DNA_BB_PROTEIN_ASA.values()),
#           (path+obj_name+"_pdiviz_dna_base.asa",                       path+obj_name+"_pdiviz_dna_base.atmasa",         stored.PDI_Viz_DNA_BASES_ASA.values()),
#           (path+obj_name+"_pdiviz_dna_base_protein.asa",               path+obj_name+"_pdiviz_dna_base_protein.atmasa", stored.PDI_Viz_DNA_BASES_PROTEIN_ASA.values()),
#           (path+obj_name+"_pdiviz_dna_wo_major_groove_protein.asa",    path+obj_name+"_pdiviz_dna_wo_major_groove_protein.atmasa", stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA.values()),
 #          (path+obj_name+"_pdiviz_dna_wo_minor_groove_protein.asa",    path+obj_name+"_pdiviz_dna_wo_minor_groove_protein.atmasa", stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA.values()),
           ]
  
  out = open(path+obj_name+"_pdiviz_stats.csv",'w')
  stats = stored.PDI_Viz_Stats[obj_name]
  out.write("OBJECT\tAREA(A^2)" + Lend)
  out.write('Complex SASA\t'                + stats['total']         + Lend +
      'Free protein SASA\t'                 + stats['protein']       + Lend +
      'Free DNA SASA\t'                     + stats['dna']           + Lend +
      'Protein backbone SASA\t'             + stats['prot_bb']       + Lend +
      'Protein side chain SASA\t'           + stats['prot_sc']       + Lend +
      'DNA backbone SASA\t'                 + stats['dna_bb']        + Lend +
      'DNA bases SASA\t'                    + stats['dna_ba']        + Lend +
      'DNA major groove SASA\t'             + stats['dna_ma']        + Lend +
      'DNA minor groove SASA\t'             + stats['dna_mi']        + Lend +
      'Buried protein surface\t'            + stats['prot-total']    + Lend +
      'Buried protein backbone surface\t'   + stats['prot-total_bb'] + Lend +
      'Buried protein side chain surface\t' + stats['prot-total_sc'] + Lend +
      'Buried DNA surface\t'                + stats['dna-total']     + Lend +
      'Buried DNA backbone surface\t'       + stats['dna-total_bb']  + Lend +
      'Buried DNA bases surface\t'          + stats['dna-total_ba']  + Lend +
      'Buried DNA major groove surface\t'   + stats['dna-total_ma']  + Lend +
      'Buried DNA minor groove surface\t'   + stats['dna-total_mi']  + Lend +
      'Buried Protein polar surface\t'      + stats['prot-total_prot_polar']  + Lend +
      'Buried Protein apolar surface\t'     + stats['prot-total_prot_apolar'] + Lend +
      'Buried DNA polar surface\t'          + stats['dna-total_dna_polar']    + Lend +
      'Buried DNA apolar surface\t'         + stats['dna-total_dna_apolar']   + Lend +
      'Interface Area\t'                    + stats['interface']     + Lend )
  out.close()

  for data in asa_l:
    out = open(data[0],'w')
    out_atm = open(data[1],'w')
    if data[-1] == 0:
      l = [(stored.PDI_Viz_OrX[item[:-2]],item[0],item[1],item[2],item[3],item[4],item[5]) for item in data[2]]
      out_atm.write('atmname\tresname\tchain\tresnum\ttotal_ASA\tbb_ASA\tsc_ASA\tmajorgroove_ASA\tminorgroove_ASA\tnogroove_ASA\t' +
                          'polar_ASA\tpolar_bb_ASA\tpolar_sc_ASA\tpolar_majorgroove_ASA\tpolar_minorgroove_ASA\tpolar_nogroove_ASA\t' +
                          'hyd_ASA\thyd_bb_ASA\thyd_sc_ASA\thyd_majorgroove_ASA\thyd_minorgroove_ASA\thyd_nogroove_ASA' + Lend)
      out.write("REMARK ID\tNAME\tRESN\tCHAIN\tRESI\tX\tY\tZ\t\tASA\tVDW"+Lend)
    else:
      l = [(stored.PDI_Viz_OrX[key],key[0],key[1],key[2],key[3],stored.PDI_Viz_T_ASA[key][4],data[2][key]) for key in list(data[2].keys())]
      out_atm.write('atmname\tresname\tchain\tresnum\ttotal_ASA\tbb_BSA\tsc_BSA\tmajorgroove_BSA\tminorgroove_BSA\tnogroove_BSA\t' +
                          'polar_BSA\tpolar_bb_BSA\tpolar_sc_BSA\tpolar_majorgroove_BSA\tpolar_minorgroove_BSA\tpolar_nogroove_BSA\t' +
                          'hyd_BSA\thyd_bb_BSA\thyd_sc_BSA\thyd_majorgroove_BSA\thyd_minorgroove_BSA\thyd_nogroove_BSA' + Lend)
      out.write("REMARK ID\tNAME\tRESN\tCHAIN\tRESI\tX\tY\tZ\t\tBSA\tVDW"+Lend)
 
    for item in sorted(l,key = lambda x: x[0]):
      xyz = stored.xyz[item[0]]
      out.write("ATOM  ") #tag
      out.write(str(item[0]).rjust(5,' ')) #id
      out.write(" ")#empty column 12
      if 'H' == item[1][0]:
        out.write(item[1].ljust(4,' '))#atom name
      else:
        out.write(" "+item[1].ljust(3,' '))#atom name
      out.write(" ") #empty altloc
      #print "|"+item[2]+"|",len(item[2])
      out.write(item[2].rjust(3,' ')) # resn
      out.write(" ") #empty column 21
      out.write(item[3]) #chain
      #resid + icode
      if (any(c.isalpha() for c in item[4])):
        out.write(item[4].rjust(5,' '))
      else:
        out.write(item[4].rjust(4,' ')+" ")

      out.write("   ")#empty columens 28-30
      #xyz coordinates
      out.write(('%.3f' % xyz[0]).rjust(8,' '))
      out.write(('%.3f' % xyz[1]).rjust(8,' '))
      out.write(('%.3f' % xyz[2]).rjust(8,' '))
      out.write(('%.3f' % item[6]).rjust(8,' '))
      out.write(('%.2f' % item[5]).rjust(6,' '))
      out.write(Lend)

      #atmasa
      total_ASA = 0.0
      bb_ASA = 0.0
      sc_ASA = 0.0
      majorgroove_ASA = 0.0
      minorgroove_ASA = 0.0
      nogroove_ASA = 0.0
      polar_ASA = 0.0
      polar_bb_ASA = 0.0
      polar_sc_ASA = 0.0
      polar_majorgroove_ASA = 0.0
      polar_minorgroove_ASA = 0.0
      polar_nogroove_ASA = 0.0
      hyd_ASA = 0.0
      hyd_bb_ASA = 0.0
      hyd_sc_ASA = 0.0
      hyd_majorgroove_ASA = 0.0
      hyd_minorgroove_ASA = 0.0
      hyd_nogroove_ASA = 0.0

      atomname = item[1]
      resname =  item[2]
      chain =    item[3]
      resnum = item[4]
      asa = item[6]

      if resname in amino_acids:
        if atomname in aa_backbone:
          found = True
          total_ASA += asa
          bb_ASA += asa
        elif atomname in aa_sidechain:
          found = True
          total_ASA += asa
          sc_ASA += asa
        if atomname in aa_polar_backbone:
          polar_ASA += asa
          polar_bb_ASA += asa
        elif atomname in aa_polar_sidechain:
          polar_ASA += asa
          polar_sc_ASA += asa
        if atomname in aa_hyd_backbone:
          hyd_ASA += asa
          hyd_bb_ASA += asa
        elif atomname in aa_hyd_sidechain:
          hyd_ASA += asa
          hyd_sc_ASA += asa

      elif resname in bases:
        resname2 = deoxybase_lookup.get(resname, resname)
        if atomname in dna_backbone:
          found = True
          total_ASA += asa
          bb_ASA += asa
        elif atomname in dna_majorgroove[resname2]:
          found = True
          total_ASA += asa
          sc_ASA += asa
          majorgroove_ASA += asa
        elif atomname in dna_minorgroove[resname2]:
          found = True
          total_ASA += asa
          sc_ASA += asa
          minorgroove_ASA += asa
        elif atomname in dna_ambiguous[resname2]:
          found = True
          total_ASA += asa
          sc_ASA += asa
          nogroove_ASA += asa
        if atomname in dna_polar_backbone:
          polar_ASA += asa
          polar_bb_ASA += asa
        elif atomname in dna_polar_majorgroove[resname2]:
          polar_ASA += asa
          polar_sc_ASA += asa
          polar_majorgroove_ASA += asa
        elif atomname in dna_polar_minorgroove[resname2]:
          polar_ASA += asa
          polar_sc_ASA += asa
          polar_minorgroove_ASA += asa
        elif atomname in dna_polar_ambigous[resname2]:
          polar_ASA += asa
          polar_sc_ASA += asa
          polar_nogroove_ASA += asa
        if atomname in dna_hyd_backbone:
          hyd_ASA += asa
          hyd_bb_ASA += asa
        elif atomname in dna_hyd_majorgroove[resname2]:
          hyd_ASA += asa
          hyd_sc_ASA += asa
          hyd_majorgroove_ASA += asa
        elif atomname in dna_hyd_minorgroove[resname2]:
          hyd_ASA += asa
          hyd_sc_ASA += asa
          hyd_minorgroove_ASA += asa
        elif atomname in dna_hyd_ambigous[resname2]:
          hyd_ASA += asa
          hyd_sc_ASA += asa
          hyd_nogroove_ASA += asa
      out_atm.write(atomname+'\t');
      out_atm.write(formatOutput((resname, chain, resnum, total_ASA, bb_ASA, sc_ASA, majorgroove_ASA, minorgroove_ASA, nogroove_ASA,
                                       polar_ASA, polar_bb_ASA, polar_sc_ASA, polar_majorgroove_ASA, polar_minorgroove_ASA, polar_nogroove_ASA,
                                       hyd_ASA, hyd_bb_ASA, hyd_sc_ASA, hyd_majorgroove_ASA, hyd_minorgroove_ASA, hyd_nogroove_ASA)) + Lend )
    out.close()
    out_atm.close()
  #print interacting residue list.
  intres = {}
  intatom = {}
  res_list = []
  not_res_list = []
  atom_list = []
  not_atom_list = []
  for pkey in stored.PDI_Viz_PROTEIN_ASA:
    reskey = pkey[1:]
    if not reskey in intres:
      #0:dbsa,1:free asa, 2:false sidechain bsa,3:false backbone bsa, 4:obj_name
      intres[reskey] = [0.0,0.0,0.0,0.0,0.0,0.0]
    intres[reskey][0] += stored.PDI_Viz_PROTEIN_m_TOTAL_d[pkey]
    intres[reskey][1] += stored.PDI_Viz_PROTEIN_ASA[pkey][5]
    intres[reskey][2] += stored.PDI_Viz_PROTEIN_BB_m_TOTAL_d[pkey]  
    intres[reskey][3] += stored.PDI_Viz_PROTEIN_BASE_m_TOTAL_d[pkey]
    intres[reskey][4] += stored.PDI_Viz_PROTEIN_NOMAJORG_DNA_m_TOTAL_d[pkey]
    intres[reskey][5] += stored.PDI_Viz_PROTEIN_NOMINORG_DNA_m_TOTAL_d[pkey]
    intatom[pkey] = stored.PDI_Viz_PROTEIN_m_TOTAL_d[pkey]
  for dkey in stored.PDI_Viz_DNA_ASA:
    reskey = dkey[1:]
    if not reskey in intres:
      #0:dbsa, 1:free asa, 2:false
      intres[reskey] = [0.0,0.0,0.0,0.0,0.0,0.0]
    intres[reskey][0] += stored.PDI_Viz_DNA_m_TOTAL_d[dkey]
    intres[reskey][1] += stored.PDI_Viz_DNA_ASA[dkey][5]
    intatom[dkey] = stored.PDI_Viz_DNA_m_TOTAL_d[dkey]

  for key in intres:
    bsa = intres[key][0]
    if bsa > 0.0:  
      free_asa = intres[key][1]
      burial_p = bsa/free_asa
      res_list.append([key[0],
                       key[1],
                       key[2],
                       bsa,
                       free_asa,
                       burial_p,
                       intres[key][2],
                       intres[key][3],
                       intres[key][4],
                       intres[key][5],
                       obj_name])
    else:
      free_asa = intres[key][1]
      not_res_list.append([key[0],key[1],key[2],bsa,free_asa,0.0, 0.0, 0.0,0.0,0.0, obj_name])
  for key in intatom:
    if intatom[key] > stored.PDI_Viz_dSASAcutoff:
      atom_list.append([key[0],key[1],key[2],key[3],intatom[key]])
    else:
      not_atom_list.append([key[0],key[1],key[2],key[3],intatom[key]])
  
  res_list.sort(key = lambda x : int(x[2]))
  res_list.sort(key = lambda x : x[1])
  
  not_res_list.sort(key = lambda x : int(x[2]))
  not_res_list.sort(key = lambda x : x[1])
  
  atom_list.sort(key = lambda x : x[0])
  atom_list.sort(key = lambda x : int(x[3]))
  atom_list.sort(key = lambda x : x[2])
  
  not_atom_list.sort(key = lambda x : x[0])
  not_atom_list.sort(key = lambda x : int(x[3]))
  not_atom_list.sort(key = lambda x : x[2])
  
  outf = open(path+obj_name+"_pdiviz_dbsa_byres.csv","w")
  outf.write("#RESN\tCHAIN\tRESI\tdBSA(A^2)\tFREE_ASA(A^2)\tBURIED_PROPORTION\tBASE_BSA\tBACKBONE_BSA\tMINORG\tMAJORG\tOBJECT"+Lend)
  for item in res_list:
    outf.write("\t".join([str(i) for i in item])+Lend)
  outf.close()
  
  outf = open(path+obj_name+"_pdiviz_noint_byres.csv","w")
  outf.write("#RESN\tCHAIN\tRESI\tdBSA(A^2)\tFREE_ASA(A^2)\tBURIED_PROPORTION\tBASE_BSA\tBACKBONE_BSA\tMAJORG\tMINORG\tOBJECT"+Lend)
  for item in not_res_list:
    outf.write("\t".join([str(i) for i in item])+Lend)
  outf.close()

  outf = open(path+obj_name+"_pdiviz_dbsa_byatom.csv","w")
  outf.write("#NAME\tRESN\tCHAIN\tRESI\tdBSA(A^2)"+Lend)
  for item in atom_list:
    outf.write("\t".join([str(i) for i in item])+Lend)
  outf.close()
  
  outf = open(path+obj_name+"_pdiviz_noint_byatom.csv","w")
  outf.write("#NAME\tRESN\tCHAIN\tRESI\tdBSA(A^2)"+Lend)
  for item in not_atom_list:
    outf.write("\t".join([str(i) for i in item])+Lend)
  outf.close()

def SaveFasta(obj_name, fname="", c80="0",retfasta=False):
  if obj_name == "":
    return
  c80 = int(c80)
  resn2one = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D',
              'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R',
              'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G',
              'PRO':'P', 'CYS':'C','HIP':'H','HID':'H',
              'A':'A','DA':'A','T':'T','DT':'T','G':'G','DG':'G','C':'C','DC':'C',
              'U':'U'}
  sel = "("+obj_name+")"
   
  it = {'lst':[],'chains': {} }
  cmd.iterate(sel,'lst.append([name,resn,resi,chain])',space=it)
  for l in it['lst']:
    if not l[3] in it['chains']:
      it['chains'][l[3]] = []
    it['chains'][l[3]].append(l[0:3])
  chains = {}
  ctype = {}
  for c in list(it['chains'].keys()):
    atoms = it['chains'][c]
    if not c in chains:
      chains[c] = {}
      hasP = False
      hasC3p = False
      hasCA= False
    for atom in atoms:
      try:
        resi = int(atom[2])
      except:
        resi = (atom[2])
      resn = resn2one[atom[1]]
      if "P" in atom[0]:
        hasP = True
      if "C3'" in atom[0]:
        hasC3p = True
      if "CA" in atom[0]:
        hasCA=True
      if not resi in chains[c]:
        chains[c][resi] = resn
    if hasP and hasC3p:
      ctype[c] = "mol:na"
    elif hasCA:
      ctype[c] = "mol:protein"
    else:
      ctype[c] = "mol:other"
  if fname == "":
    fname = obj_name+"_pdiviz.fasta"
  if not retfasta: 
    f = open(fname,"w")
    for c in sorted(chains.keys()):
      chain = chains[c]
      fasta = []
      i=0
      for resi in sorted(chain.keys()):
        fasta.append(chain[resi])
        i+=1
        if ((c80 == 1) and ((i % 80) == 0)):
          fasta.append("\n")
      fasta = "".join(fasta)
      f.write(">"+obj_name+"_"+c+" "+ctype[c]+"\n")
      f.write(fasta+"\n")
    f.close()
  else:
    fastastr = []
    for c in sorted(chains.keys()):
      chain = chains[c]
      fasta = []
      i=0
      for resi in sorted(chain.keys()):
        fasta.append(chain[resi])
      fasta = "".join(fasta)
      fastastr.append([obj_name,c,ctype[c],fasta])
    return fastastr
cmd.extend("PDI_Viz_SaveFasta",SaveFasta)
stored.PDI_Viz_SaveFasta = SaveFasta

################################################################################
#holds the IDs of each selection for later use
stored.S_IDs = {}


def GenSeleFromResnNameDict(obj_name,resn_name_dict):
  templist = []
  for resn in resn_name_dict:
    l = resn_name_dict[resn]
    if len(l) > 0:
      templist.append("(resn "+resn+" and name "+','.join(l)+") ")
  return "object "+obj_name + ' and ('+' or '.join(templist)+')'

def SeparateComplex(obj_name):
  """
  Separates a single protein-dna complex into protein, dna, protein backbone and sidechain, and dna backbone and bases selections.
  Also separates in major/minor groove atoms.
  Also separates base DNA atoms as h bond donors, acceptors and thymine methyls.
  Now also identifies polar and apolar surfaces in protein and DNA
  obj_name = object to separate (string)
  """

  p_bb = obj_name + ' and name '+','.join(stored.PDI_Viz_prot_backbone)
  p_sc = obj_name + ' and name '+','.join(stored.PDI_Viz_prot_base)
  dna_bb = obj_name + ' and name '+','.join(stored.PDI_Viz_dna_backbone)
  dna_ba = obj_name + ' and name '+','.join(stored.PDI_Viz_dna_base)

  dna_ma = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_dna_majorgroove)
  dna_mi = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_dna_minorgroove)
  dna_amb = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_dna_ambi)

  dna_polar_ba = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_dna_polar_ba)
  dna_polar_bb = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_dna_polar_bb)
  dna_polar = '(('+dna_polar_ba+') or ('+dna_polar_bb+'))'
  dna_apolar_ba = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_dna_apolar_ba)
  dna_apolar_bb = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_dna_apolar_bb)
  dna_apolar = '(('+dna_apolar_ba+') or ('+dna_apolar_bb+'))'
  
  prot_polar = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_prot_polar)
  prot_apolar = GenSeleFromResnNameDict(obj_name,stored.PDI_Viz_prot_apolar)
  
  
  #holds id keys for atom selections
  idkey = {'p_bb': [],
           'p_sc': [],
           'p_polar':[],
           'p_apolar':[],
           'dna_bb': [],
           'dna_ba': [],
           'dna_ma':[],
           'dna_mi':[],
           'dna_amb':[],
           'dna_polar':[],
           'dna_apolar':[]}

  #creates selection with protein and DNA separated
  cmd.select(stored.PDI_Viz_p_bb, p_bb )
  cmd.select(stored.PDI_Viz_p_sc, p_sc )
  
  cmd.select(stored.PDI_Viz_DNA_bb, dna_bb )
  cmd.select(stored.PDI_Viz_DNA_base, dna_ba )
  
  cmd.select(stored.PDI_Viz_Prot, obj_name + ' and bychain name CA')
  cmd.select(stored.PDI_Viz_DNA, obj_name + ' and bychain name P')
  
  cmd.select(stored.PDI_Viz_DNA_majorg, dna_ma)
  cmd.select(stored.PDI_Viz_DNA_minorg, dna_mi)
  cmd.select(stored.PDI_Viz_DNA_amb, dna_amb)
  
  cmd.select(stored.PDI_Viz_DNA_polar,dna_polar)
  cmd.select(stored.PDI_Viz_DNA_apolar,dna_apolar)
  
  cmd.select(stored.PDI_Viz_p_polar, prot_polar )
  cmd.select(stored.PDI_Viz_p_apolar, prot_apolar )

  #stores the atom id key for each selection
  cmd.iterate(stored.PDI_Viz_p_bb, 'p_bb.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_p_sc, 'p_sc.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_DNA_bb, 'dna_bb.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_DNA_base, 'dna_ba.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_DNA_majorg, 'dna_ma.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_DNA_minorg, 'dna_mi.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_DNA_amb, 'dna_amb.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_DNA_polar, 'dna_polar.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_DNA_apolar, 'dna_apolar.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_p_polar, 'p_polar.append('+id_t+')',space=idkey)
  cmd.iterate(stored.PDI_Viz_p_apolar, 'p_apolar.append('+id_t+')',space=idkey)
  stored.S_IDs[obj_name] = idkey

def DeleteComplex():
  """
  Deletes selections created by SeparateComplex
  No arguments needed.
  """
  cmd.delete(stored.PDI_Viz_p_bb)
  cmd.delete(stored.PDI_Viz_p_sc)
  cmd.delete(stored.PDI_Viz_DNA_bb)
  cmd.delete(stored.PDI_Viz_DNA_base)
  cmd.delete(stored.PDI_Viz_Prot)
  cmd.delete(stored.PDI_Viz_DNA)
  cmd.delete(stored.PDI_Viz_DNA_majorg)
  cmd.delete(stored.PDI_Viz_DNA_minorg)
  cmd.delete(stored.PDI_Viz_DNA_amb)
  cmd.delete(stored.PDI_Viz_DNA_polar)
  cmd.delete(stored.PDI_Viz_DNA_apolar)
  cmd.delete(stored.PDI_Viz_p_polar)
  cmd.delete(stored.PDI_Viz_p_apolar)

  cmd.delete(stored.PDI_Viz_bGTZeroProt)
  cmd.delete(stored.PDI_Viz_qGTZeroProt)
  cmd.delete(stored.PDI_Viz_bGTZeroDNA)
  cmd.delete(stored.PDI_Viz_qGTZeroDNA)
  cmd.delete(stored.PDI_Viz_bGTZeroProt2)
  cmd.delete(stored.PDI_Viz_qGTZeroProt2)
  cmd.delete(stored.PDI_Viz_bGTZeroDNA2)
  cmd.delete(stored.PDI_Viz_qGTZeroDNA2)

cmd.extend("PDI_Viz_SeparateComplex",SeparateComplex)
cmd.extend("PDI_Viz_DeleteComplex",DeleteComplex)
stored.PDI_Viz_SeparateComplex = SeparateComplex
stored.PDI_Viz_DeleteComplex = DeleteComplex
################################################################################
def fixname(name):
  m = re.search("O[1-5]P",name)
  try:
    if m.group(0) == name:
      name = "OP"+name[1]
      return name
  except:
    pass
  if name == 'C5M':
      return 'C7'
  else:
    if "*" in name:
      return name.replace('*','\'')
  return name
def FixName(obj_name=""):
  """
  Fixes older atom names with *.
  obj_name = (string) (default "")
  """
  f = {"Fname" : fixname}
  cmd.alter(obj_name,"name=Fname(name)",space=f)

def FixVDW(obj_name = ""):
  """
  Fixes VDW radii of the given object or all structures.
  obj_name = (string) (default "")
  """
  if obj_name != "":
    obj_name = obj_name + " and "
  res_dict = stored.PDI_Viz_res_dict
  for resn in list(res_dict.keys()):
    names = list(res_dict[resn].keys())
    for name in names:
      #print obj_name+'resn "'+resn+'" and name "'+name+'"'
      cmd.alter(obj_name+'resn "'+resn+'" and name "'+name+'"','vdw = stored.PDI_Viz_res_dict["'+resn+'"]["'+name+'"]')
  cmd.rebuild()
  print("Atom radii fixed.")

def RemoveAltAtoms(obj_name = ""):
  """
  Removes the lowest occupancy atoms from the current structure or from the object name give.
  Assumes that the PDB is consistent in the alternate location style, with alternate positions in the same
  chain or separated in consecutive chains.
  obj_name = (string) (default "")
  """
  if obj_name != "":
    obj_name = obj_name + " and "

  stored.PDI_Viz_AltAtoms = []
  cmd.iterate(obj_name+' not alt ""','stored.PDI_Viz_AltAtoms.append((ID,name,resi,chain,alt,q))')

  chain_dict = {}
  count = 0
  alt_a = stored.PDI_Viz_AltAtoms
  for atom in alt_a:
    name  = atom[1]
    resi  = atom[2]
    chain = atom[3]
    if( not chain in chain_dict):
      chain_dict[chain] = {}
    if( not resi in chain_dict[chain]):
      chain_dict[chain][resi] = {}
    if (not name in chain_dict[chain][resi]):
      chain_dict[chain][resi][name] = []

    chain_dict[chain][resi][name].append(atom)
  AltInChain = True
  ValidAlternateFormat = True
  for chain in list(chain_dict.values()):
    if(AltInChain):
      for res in list(chain.values()):
        if(AltInChain):
          for alt_atoms in list(res.values()):
            if (len(alt_atoms) < 2):
              AltInChain = False
              break
            alt_atoms.sort(key = lambda x: x[4], reverse = True)
            alt_atoms.sort(key = lambda x: x[5])
            for atom in alt_atoms[:-1]:
              print("Removing atom: "+obj_name+'id '+str(atom[0]))
              cmd.remove(obj_name+'id '+str(atom[0]))
              count+=1

        else:
          break
    else:
      break

  if(not AltInChain):
    #assumes chains are in order
    #can only deal with 2 altloc chains
    alt_a.sort(key = lambda x: x[0])
    chain_dict = {}
    for atom in alt_a:
      chain = atom[3]
      if( not chain in chain_dict):
        chain_dict[chain] = []
      chain_dict[chain].append(atom)
    chains = list(chain_dict.keys())
    chain_avg = []
    for key in chains:
      t = [item[5] for item in chain_dict[key]]
      chain_avg.append((key,round(sum(t)/len(chain_dict[key]),3)))
    eq_chains = []
    i = 0
    j = 1
    l = len(chains)
    while j<l:
      eq_chains.append([chain_avg[i],chain_avg[j]])
      i+=2
      j+=2
    for pair in eq_chains:
      pair.sort(key = lambda x: x[0],reverse = True)
      pair.sort(key = lambda x: x[1])
      count+= len(chain_dict[pair[0][0]])
      print(pair[0][0])
      cmd.remove(obj_name+' chain '+pair[0][0])

  if(ValidAlternateFormat):
    cmd.alter('all','alt=""')
  print(count,"alternate atoms removed")
  return not AltInChain

cmd.extend ("PDI_Viz_FixVDW", FixVDW)
cmd.extend ("PDI_Viz_RemoveAltAtoms", RemoveAltAtoms)
stored.PDI_Viz_FixVDW = FixVDW
stored.PDI_Viz_RemoveAltAtoms = RemoveAltAtoms

################################################################################
def ReloadASACalc(obj_name):
  if(obj_name in stored.PDI_Viz_Results):
      r = stored.PDI_Viz_Results[obj_name]
      stored.PDI_Viz_T_ASA = r[0]
      stored.PDI_Viz_PROTEIN_ASA = r[1]
      stored.PDI_Viz_DNA_ASA = r[2]
      stored.PDI_Viz_DNA_BB_ASA = r[3]
      stored.PDI_Viz_DNA_BB_PROTEIN_ASA = r[4]
      stored.PDI_Viz_DNA_BASES_ASA = r[5]
      stored.PDI_Viz_DNA_BASES_PROTEIN_ASA = r[6]
      stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA = r[7]
      stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA = r[8]
      stored.PDI_Viz_DNA_MAJOR_GROOVE_ASA = r[9]
      stored.PDI_Viz_DNA_MINOR_GROOVE_ASA = r[10]
      return True
  else:
    return False

def PrintOutputFiles(obj_name = "",path = './',savesession = "-1",fasta="0"):
  """
  Prints *.asa, *.atmasa, dna color and stats files.
  obj_name = Name of object to analyze (string)
  path = Path to store output files (string) (default "./")
  prot = Print protein aminoacid color (string) (0 or 1, default 0)
  dna_byres = Print dna colors (int) (0 or 1, default 1)
  dna_byatom = Prints dna colors by atom (int) (0 or 1, default 0)
  savesession = Save visualization structure session. -1 Disables save session,
                0 saves all structures,  1-9 Saves the selected visualization. (int) (-1 to 9, default -1)
  """

  if(obj_name == ""):
    #assumes user only loaded a single pdb
    if(len(cmd.get_object_list()) > 0):
      obj_name = cmd.get_object_list()[0]
    else:
      print("Nothing loaded.")

  if not IsObjProtDNAcomplex(obj_name):
    print("Structure "+obj_name+" is not a Protein/DNA complex.")
    return
  if IsObjMultiState(obj_name):
    print("Multi state objects currently unsupported. Please split the model.")
    return
  cmd.wizard("message","Saving data for "+obj_name)
  print("Begin file output for "+obj_name)
  cmd.set('suspend_updates', 'on')
  SetUpEnv(obj_name)
  view = cmd.get_view()
  SeparateComplex(obj_name)
  StoreOriginalAtomID(obj_name)
  if(not ReloadASACalc(obj_name)):
    CalculateSplitASA(obj_name)
  else:
    print("Precalculated data found.")
  CalculateDelta()
  StoreAreas(obj_name)
  WriteSplitASA(obj_name,path)
  #ColorOut(obj_name,path,int(prot),int(dna_byres),int(dna_byatom))
  DeleteComplex()
  cmd.set_view(view)
  cmd.set('suspend_updates', 'off')
  cmd.wizard()
  if int(savesession) != -1:
    Visualize(obj_name,mode=savesession,recalc='0')
    cmd.save(path+obj_name+"_pdiviz_"+savesession+".pse", obj_name+"*")
    cmd.delete(obj_name+"_*" )
  if int(fasta) != 0:
    SaveFasta(obj_name,path+obj_name+".fasta","1")
  print("DATA SAVED!")

cmd.extend('PDI_Viz_PrintFiles',PrintOutputFiles)
stored.PDI_Viz_PrintFiles = PrintOutputFiles

def Visualize(obj_name = "",mode = "0" ,bg=sd_set[0], dna_c=sd_set[1], dna_bb_c=sd_set[2],dna_ba_c=sd_set[3],dna_both=sd_set[4],T=sd_set[5],grad0 = sd_set[6],grad1= sd_set[7],grad2=sd_set[8],rohs_color=sd_set[9],recalc='1'):
  """
  Creates a dna/protein interaction visualization of the current or selected object.
  mode = Selects drawing mode (string) (1 or 2, default 2)
  obj_name = object to use (string)
  bg = background color rgb tuple (string)
  dna_c = DNA initial color RGB tuple (string)
  dna_bb_c = DNA backbone/protein interaction color rgb tuple (string)
  dna_ba_c = DNA base/protein interaction color rgb tuple (string)
  dna_both = DNA bb and base / protein interaction color rgb tuple (string)
  T = surface gransparency color rgb tuple (string)
  grad0, grad1, grad2 = strings containing tuples of 9 float values (0.0/1.0) for low,mid,high ramp
  dbg = Enable file output (string) (0 or 1, default 0)
  dbg_path = File output folder (string) (default "./")
  out_color = Enable DNA color output. (string) (0 or 1, default 0)
  """
  view = cmd.get_view()

  cmd.set_color("PDI_Viz_bg", [float(item) for item in bg[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_c", [float(item) for item in dna_c[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_bb_c", [float(item) for item in dna_bb_c[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_ba_c", [float(item) for item in dna_ba_c[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_both", [float(item) for item in dna_both[1:-1].split(',')])

  cmd.set_color("PDI_Viz_dna_noint",[float(item) for item in rohs_color[1:-1].split(',')[0:3]])
  cmd.set_color("PDI_Viz_dna_acceptor",[float(item) for item in rohs_color[1:-1].split(',')[3:6]])
  cmd.set_color("PDI_Viz_dna_donor",[float(item) for item in rohs_color[1:-1].split(',')[6:9]])
  cmd.set_color("PDI_Viz_dna_thymine",[float(item) for item in rohs_color[1:-1].split(',')[9:12]])
  cmd.set_color("PDI_Viz_dna_neutral",[float(item) for item in rohs_color[1:-1].split(',')[12:15]])

  #set interaction sequence colors
  cmd.set_color('dna_DA', (1.00,0.00,1.00))
  cmd.set_color('dna_DH', (0.50,0.50,1.00))
  cmd.set_color('dna_AH', (1.00,0.50,0.50))
  cmd.set_color('dna_DAH',(0.83,0.71,0.83))
  #thymine exclusive
  cmd.set_color('dna_DT', (0.43,0.64,0.66))
  cmd.set_color('dna_AT', (1.00,0.02,0.32))
  cmd.set_color('dna_HT', (0.98,1.00,0.35))
  cmd.set_color('dna_DAT',(0.11,1.00,0.56))
  cmd.set_color('dna_DTH',(0.73,0.06,0.93))
  cmd.set_color('dna_ATH',(1.00,0.55,0.13))
  cmd.set_color('dna_DATH',(0.74,1.0,0.85))

  T = float(T)
  prot_grad = (grad0,grad1,grad2)
  if(obj_name == ""):
    #assumes user only loaded a single pdb
    if(len(cmd.get_object_list()) > 0):
      obj_name = cmd.get_object_list()[0]
    else:
      print("Nothing loaded.")
      return False

  if not IsObjProtDNAcomplex(obj_name):
    print("Structure "+obj_name+" is not a Protein/DNA complex.")
    return False
  if IsObjMultiState(obj_name):
    print("Multi state objects currently unsupported. Please split the model")
    return False

  cmd.delete(obj_name+'_*')
  SetUpEnv(obj_name)
  SeparateComplex(obj_name)
  StoreOriginalAtomID(obj_name)

  if(recalc == '0'):
    if(not ReloadASACalc(obj_name)):
      DeleteComplex()
      return False
  else:
    CalculateSplitASA(obj_name)
  CalculateDelta()
  StoreAreas(obj_name)
  Mode(mode,obj_name,
       "PDI_Viz_bg","PDI_Viz_dna_c","PDI_Viz_dna_bb_c","PDI_Viz_dna_ba_c","PDI_Viz_dna_both",
       T,prot_grad,
       "PDI_Viz_dna_noint","PDI_Viz_dna_acceptor","PDI_Viz_dna_donor","PDI_Viz_dna_thymine","PDI_Viz_dna_neutral")
  Show(obj_name)
  DeleteComplex()

  cmd.set_view(view)
  return True

def Visualize2(obj_name = "",mode = "2" ,bg=sd_set[0], dna_c=sd_set[1], dna_bb_c=sd_set[2],dna_ba_c=sd_set[3],dna_both=sd_set[4],T=sd_set[5],grad0 = sd_set[6],grad1= sd_set[7],grad2=sd_set[8],rohs_color=sd_set[9],recalc='1',dbg='0',dbg_path="",out_color='0'):
  """
  Creates a dna/protein interaction visualization of the current or selected object.
  mode = Selects drawing mode (string) (1 or 2, default 2)
  obj_name = object to use (string)
  bg = background color rgb tuple (string)
  dna_c = DNA initial color rgb tuple (string)
  dna_bb_c = DNA backbone/protein interaction color rgb tuple (string)
  dna_ba_c = DNA base/protein interaction color rgb tuple (string)
  dna_both = DNA bb and base / protein interaction color rgb tuple (string)
  T = surface gransparency color rgb tuple (string)
  grad0, grad1, grad2 = strings containing tuples of 9 float values (0.0/1.0) for low,mid,high ramps
  dbg = Enable file output (string) (0 or 1, default 0). Deprecated.
  dbg_path = File output folder (string) (default "./") Deprecated.
  out_color = Enable DNA color output. (string) (0 or 1, default 0) Deprecated.
  """
  view = cmd.get_view()

  cmd.set_color("PDI_Viz_bg", [float(item) for item in bg[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_c", [float(item) for item in dna_c[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_bb_c", [float(item) for item in dna_bb_c[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_ba_c", [float(item) for item in dna_ba_c[1:-1].split(',')])
  cmd.set_color("PDI_Viz_dna_both", [float(item) for item in dna_both[1:-1].split(',')])

  cmd.set_color("PDI_Viz_dna_noint",[float(item) for item in rohs_color[1:-1].split(',')[0:3]])
  cmd.set_color("PDI_Viz_dna_acceptor",[float(item) for item in rohs_color[1:-1].split(',')[3:6]])
  cmd.set_color("PDI_Viz_dna_donor",[float(item) for item in rohs_color[1:-1].split(',')[6:9]])
  cmd.set_color("PDI_Viz_dna_thymine",[float(item) for item in rohs_color[1:-1].split(',')[9:12]])
  cmd.set_color("PDI_Viz_dna_neutral",[float(item) for item in rohs_color[1:-1].split(',')[12:15]])

  T = float(T)
  prot_grad = (grad0,grad1,grad2)

  if not IsObjProtDNAcomplex(obj_name):
    print("Structure "+obj_name+" is not a Protein/DNA complex.")
    return
  if IsObjMultiState(obj_name):
    print("Multi state objects currently unsupported. Please split the model")
    return

  if(recalc == '0'):
    if(obj_name in stored.PDI_Viz_Results):
      r = stored.PDI_Viz_Results[obj_name]
      stored.PDI_Viz_T_ASA = r[0]
      stored.PDI_Viz_PROTEIN_ASA = r[1]
      stored.PDI_Viz_DNA_ASA = r[2]
      stored.PDI_Viz_DNA_BB_ASA = r[3]
      stored.PDI_Viz_DNA_BB_PROTEIN_ASA = r[4]
      stored.PDI_Viz_DNA_BASES_ASA = r[5]
      stored.PDI_Viz_DNA_BASES_PROTEIN_ASA = r[6]
      stored.PDI_Viz_PROT_DNA_NO_MAJOR_GROOVE_ASA = r[7]
      stored.PDI_Viz_PROT_DNA_NO_MINOR_GROOVE_ASA = r[8]
      stored.PDI_Viz_DNA_MAJOR_GROOVE_ASA = r[9]
      stored.PDI_Viz_DNA_MINOR_GROOVE_ASA = r[10]
      Colorize(mode,obj_name,
       "PDI_Viz_bg","PDI_Viz_dna_c","PDI_Viz_dna_bb_c","PDI_Viz_dna_ba_c","PDI_Viz_dna_both",
       T,prot_grad,
       "PDI_Viz_dna_noint","PDI_Viz_dna_acceptor","PDI_Viz_dna_donor","PDI_Viz_dna_thymine","PDI_Viz_dna_neutral")
  cmd.set_view(view)

def IsObjMultiState(obj_name):
  return (cmd.count_states(obj_name) > 1)

def IsObjProtDNAcomplex(obj_name):
  for item in ["_prot_obj_m2","_dna_int_m2","_dna_obj_m2","_prot_obj_g3","_dna_obj_g3"]:
    if item in obj_name:
      return False
  #if obj_name == "":
    #return False
  #print 1,"|"+obj_name+"|"
  p_a = cmd.select('prot_test',"object "+ obj_name+' and bychain name CA')
  #print 2,"|"+obj_name+"|"
  cmd.delete('prot_test')
  #print 3,"|"+obj_name+"|"
  d_a = cmd.select('dna_test',"object "+ obj_name+' and bychain name P')
  #print 4,"|"+obj_name+"|"
  cmd.delete('dna_test')
  #print 5,"|"+obj_name+"|"
  return (p_a > 0 and d_a > 0)

cmd.extend("PDI_Viz_Visualize",Visualize)
stored.PDI_Viz_Visualize = Visualize

################################################################################

################################################################################
'''
http://pymolwiki.org/index.php/cgo_arrow

(c) 2013 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
       color='blue red', name=''):
  '''
DESCRIPTION

  Create a CGO arrow between two picked atoms.

ARGUMENTS

  atom1 = string: single atom selection or list of 3 floats {default: pk1}

  atom2 = string: single atom selection or list of 3 floats {default: pk2}

  radius = float: arrow radius {default: 0.5}

  gap = float: gap between arrow tips and the two atoms {default: 0.0}

  hlength = float: length of head

  hradius = float: radius of head

  color = string: one or two color names {default: blue red}

  name = string: name of CGO object
  '''
  from chempy import cpv

  radius, gap = float(radius), float(gap)
  hlength, hradius = float(hlength), float(hradius)

  try:
    color1, color2 = color.split()
  except:
    color1 = color2 = color
  color1 = list(cmd.get_color_tuple(color1))
  color2 = list(cmd.get_color_tuple(color2))

  def get_coord(v):
    if not isinstance(v, str):
      return v
    if v.startswith('['):
      return cmd.safe_list_eval(v)
    return cmd.get_atom_coords(v)

  xyz1 = get_coord(atom1)
  xyz2 = get_coord(atom2)
  normal = cpv.normalize(cpv.sub(xyz1, xyz2))

  if hlength < 0:
    hlength = radius * 3.0
  if hradius < 0:
    hradius = hlength * 0.6

  if gap:
    diff = cpv.scale(normal, gap)
    xyz1 = cpv.sub(xyz1, diff)
    xyz2 = cpv.add(xyz2, diff)

  xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

  obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
     [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
     [1.0, 0.0]
  #obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
     #[cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
     #[1.0, 0.0]

  if not name:
    name = cmd.get_unused_name('arrow')

  cmd.load_cgo(obj, name)
stored.cgo_arrow = cgo_arrow
cmd.extend('DRSASA_cgo_arrow', cgo_arrow)
################################################################################
def ReadDRSASAIntTable(fname):
  lines = open(fname).readlines()
  if("<-" in lines[0]):
    direction = "col"
  elif ("->" in lines[0]):
    direction = "row"
  else:
    raise "incompatible file?"
  tag = lines[0].split("\t")[0]
  columns = lines[0].strip().split("\t")[1:]
  sumv = []
  #print columns
  columns = [ item.split("/") for item in columns]
  rows = []
  matrix = []
  for line in lines[1:]:
    line = line.strip()
    tokens = [item.strip() for item in line.split("\t")]
    if len(tokens) < 3:
      continue
    rlbl = tokens[0].split("/")
    #print rlbl
    rows.append(rlbl)
    matrix.append([float(value) for value in tokens[1:]])
  if direction == "row":
    for row in range(len(matrix)):
      sumv.append(sum(matrix[row]))
  elif direction == "col":
    sumv = [0.0 for i in range(len(matrix[0]))]
    for row in range(len(matrix)):
      for col in range(len(matrix[row])):
        sumv[col]+=matrix[row][col]
    
  return columns,rows,matrix,direction,sumv,tag

def GenSeleResn(resn,chain,resi):
  if resn in ["DG","DA","A","G"]:
    atom = "(name N9) and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  elif resn in ["DC","DT","C","T","U"]:
    atom = "(name N1) and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  else:
    atom = "(name CA) and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  
  resn = "(resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  return atom,resn

def GenSeleAtom(name,resn,chain,resi):
  atom = "(name "+name+") and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  return atom

def DrawDRSASA_ByRes(fname):
  cols,rows,matrix,direc,sumv,tag = ReadDRSASAIntTable(fname)
  mtrx = np.array(matrix).flatten()
  mtrx = mtrx[mtrx != 0]
  mvalue = np.max(mtrx)
  minvalue = np.min(mtrx)
  mv075 = mvalue*0.75
  median = np.median(mtrx)
  avg = np.average(mtrx)
  cutoff = mvalue*0.05
  vsumnp = np.array(sumv)
  vsumnp = vsumnp[vsumnp != 0]
  maxvsum = np.max(vsumnp)
  minvsum = np.min(vsumnp)
  print(median,avg,cutoff,minvalue)
  Rl = []
  Gl = []
  Bl = []
  cmd.set('suspend_updates', 'on')
  view = cmd.get_view()
  for row in range(len(matrix)):
    for col in range(len(matrix[row])):
      v = float(matrix[row][col])
      if v < cutoff:
        continue
      R = 1.8*v/mvalue
      
      rid = [str(item) for item in rows[row]]
      cid = [str(item) for item in cols[col]]
      sele1,resn1 = GenSeleResn(rid[0],rid[1],rid[2])
      sele2,resn2 = GenSeleResn(cid[0],cid[1],cid[2])
      print((sele1,sele2))
      cmd.select("atom1",sele1)
      cmd.select("atom2",sele2)
      d = cmd.get_distance("atom1","atom2")
      #name = "/".join(rid)+"_"+"/".join(cid)
      name = "a_"+str(row)+"_"+str(col)
      if direc == "row":
        resn_s = resn2
        resn_o = resn1
        vsum = sumv[row]
      elif direc == "col":
        resn_s = resn1
        resn_o = resn2
        vsum = sumv[col]
      vsumR = vsum/maxvsum
      if (v >= mv075):
        color = "blue blue"
        Rl.append(resn_s)
        cmd.color("blue",resn_s)
      elif ( v > median):
        color = "green green"
        Gl.append(resn_s)
        if not resn_s in Rl:
          cmd.color("green",resn_s)
      else:
        color = "red red"
        if (not resn_s in Rl) and (not resn_s in Gl):
          cmd.color("red",resn_s)
      coloro = "drsasa_o_"+str(row)+"_"+str(col)
      cmd.set_color(coloro,(int(255*vsumR),0, 255 - int(255*vsumR)))
      cmd.color(coloro,resn_o)
      if direc == "row":
        cgo_arrow("atom1","atom2",0.1, d*0.25 , d*0.25 , R, color, name)
      elif direc == "col":
        cgo_arrow("atom2","atom1",0.1, d*0.25 , d*0.25 , R, color, name)
  cmd.set_view(view)
  cmd.set('suspend_updates', 'off')
cmd.extend("DRSASA_DrawInter_ByRes",DrawDRSASA_ByRes)

def DrawDRSASA_ByAtom(fname,mode='0'):
  mode = int(mode)
  cols,rows,matrix,direc,sumv,tag = ReadDRSASAIntTable(fname)
  mtrx = np.array(matrix).flatten()
  mtrx = mtrx[mtrx != 0]
  mvalue = np.max(mtrx)
  minvalue = np.min(mtrx)
  mv075 = mvalue*0.75
  median = np.median(mtrx)
  avg = np.average(mtrx)
  cutoff = mvalue*0.05
  vsumnp = np.array(sumv)
  vsumnp = vsumnp[vsumnp != 0]
  maxvsum = np.max(vsumnp)
  minvsum = np.min(vsumnp)
  
  
  print(median,avg,cutoff,minvalue)
  Rl = []
  Gl = []
  Bl = []
  
  view = cmd.get_view()
  cmd.set('suspend_updates', 'on')
  cmd.set("sphere_scale",0.3)
  for row in range(len(matrix)):
    for col in range(len(matrix[row])):
      v = float(matrix[row][col])
      if v < cutoff:
        continue
      R = 1.8*v/mvalue
      nR = v/mvalue
      colorn = "drsasa_"+str(row)+"_"+str(col)
      coloro = "drsasa_o_"+str(row)+"_"+str(col)
      cmd.set_color(colorn,(int(255*nR),0, 255 - int(255*nR)))
      rid = [str(item) for item in rows[row]]
      cid = [str(item) for item in cols[col]]
      print((rid,cid))
      sele1 = GenSeleAtom(rid[0],rid[1],rid[2],rid[3])
      sele2 = GenSeleAtom(cid[0],cid[1],cid[2],cid[3])
      print((sele1,sele2))
      cmd.select("atom1",sele1)
      cmd.select("atom2",sele2)
      d = cmd.get_distance("atom1","atom2")
      if direc == "row":
        atom_s = sele2
        atom_o = sele1
        vsum = sumv[row]
      elif direc == "col":
        atom_s = sele1
        atom_o = sele2
        vsum = sumv[col]
      vsumR = vsum/maxvsum
      cmd.set_color(coloro,(int(255*vsumR),0, 255 - int(255*vsumR)))
      #name = "/".join(rid)+"_"+"/".join(cid)
      name = "a_"+str(row)+"_"+str(col)
      cmd.show("sphere",atom_s)
      cmd.show("sphere",atom_o)
      if mode == 0:
        if (v >= mv075):
          color = "red red"
          Rl.append(atom_s)
          cmd.color("red",atom_s)
        elif ( v > median):
          color = "green green"
          Gl.append(atom_s)
          if not atom_s in Rl:
            cmd.color("green",atom_s)
        else:
          color = "blue blue"
          if (not atom_s in Rl) and (not atom_s in Gl):
            cmd.color("blue",atom_s)
      elif mode == 1:
        color = colorn+" "+colorn
        cmd.color(colorn,atom_s)
      cmd.color(coloro,atom_o)
      if direc == "row":
        if mode == 0:
          cgo_arrow("atom1","atom2",0.1, d*0.25 , d*0.25 , R, color, name)
        elif mode == 1:
          cgo_arrow("atom1","atom2",0.1, d*0.25 , d*0.25 , 0.2, color, name)
      elif direc == "col":
        if mode == 0:
          cgo_arrow("atom2","atom1",0.1, d*0.25 , d*0.25 , R, color, name)
        elif mode == 1:
          cgo_arrow("atom2","atom1",0.1, d*0.25 , d*0.25 , 0.2, color, name)
  cmd.set_view(view)
  cmd.set('suspend_updates', 'off')
cmd.extend("DRSASA_DrawInter_ByAtom",DrawDRSASA_ByAtom)

def DR_SASA_ColorSASAE(objname,efile):
  ftokens = [item.strip().split("\t") for item in open(efile).readlines()]
  stored.repdict = {}
  elist = []
  for tokens in ftokens:
    stored.repdict[(tokens[1],tokens[2],tokens[4],tokens[3])] = float(tokens[5])
    elist.append(float(tokens[5]))
  cmd.alter(objname,"b = stored.repdict["+id_t+"]")
  median = elist[len(elist)/2]
  maxv = max(elist)
  minv = min(elist)
  color_b(objname,gradient='bwr')
  print(("Values: ",str(minv)+" "+str(median)+" "+str(maxv)))
cmd.extend("DRSASA_ColorESASA",DR_SASA_ColorSASAE)
################################################################################
#color_b source code

# rlc color_b.py version 8.0
# Copyright (c) 2004 Robert L. Campbell
# added user defined colors for 3 color ramp- Mark A. Wall
# added selection string to color names to avoid overlap between different selections
#
# Modified by Andreas Schueller <aschueller@bio.puc.cl>
# MOD_VERSION: 0.1 [2013-06-04]  Added item='bq' functionality


# main function called from within PyMOL

def color_b(selection='all',item='b',mode='hist',gradient='bgr',nbins=11,sat=1.,value=1.,minimum='',maximum='',user_rgb='',debug=0):
  """

  AUTHOR

    Robert L. Campbell with enhancements from James Stroud
    Modified by Andreas Schueller <aschueller@bio.puc.cl>
    MOD_VERSION: 0.1 [2013-06-04]

  USAGE

    color_b selection='sel',item='b', 'q', 'partial_charge' or 'formal_charge'
      gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' or
      'rw' or 'wr' or 'ry' or 'yr' or 'gw' or 'wg' or 'bw' or wb' or
      'gy' or 'yg' or 'gray' or 'reversegray' or 'user'

      mode='hist' or 'ramp' (default is 'hist')

      [minimum=''],[maximum=20.],

      nbins=11, sat=1.0, value=1.0,

      user_rgb = '(r1,g1,b1,r2,g2,b2,r3,g3,b3') [for use with gradient=user]

      The "item" argument allows specifying 'b', 'q', 'index', 'partial_charge'
      or 'formal_charge'as the item to color on.  The "color_q" function
      is really just the same as "color_b item=q".  Using item=index is
      similar to using the built-in "spectrum" command.

      This function allows coloring of a selection as a function of
      B-value or occupancy, following a gradient of colours.  The
      gradients can be:

      'bgr': blue -> green   -> red
      'rgb': red  -> green   -> blue
      'bwr': blue -> white   -> red
      'rwb': red  -> white   -> blue
      'bmr': blue -> magenta -> red
      'rmb': red  -> magenta -> blue
      'rw' : red -> white
      'wr' : white -> red
      'ry' : red -> yellow
      'yr' : yellow -> red
      'gw' : green -> white
      'wg' : white -> green
      'bw' : blue -> white
      'wb' : white -> blue
      'gy' : green -> yellow
      'yg' : yellow -> green
      'gray' : black -> white
      'reversegray' : white -> black
      'user' : user defined in this script

      ('rainbow' and 'reverserainbow' can be used as synonyms for
      'bgr' and 'rgb' respectively and 'grey' can be used as a synonym for 'gray').

      User-defined gradients are entered on the command line in
      parentheses as either integers between 0 and 255 or floats between
      0 and 1.  If any one value is larger than 1, then it is assumed
      that all are being entered as integers between 0 and 255.  Hence one can type:

      color_b selection, gradient=user, user_rgb=(0,0,1, 0,.5,1., 1.,.5,0.)

        or

      color_b selection, gradient=user, user_rgb=(0,0,255, 0,128,255, 255,128,0.)

      The division of B-value ranges can be in either of two modes: 'hist' or
      'ramp'. 'hist' is like a histogram (equal-sized B-value increments
      leading to unequal numbers of atoms in each bin). 'ramp' as a ramp
      of B-value ranges with the ranges chosen to provide an equal number
      of atoms in each group.

      You can also specify the lower or upper limits of the data used to determine
      the color bins (minimum,maximum). e.g. color_b my_molecule, minimum=15., maximum=25.

      You can also specify the saturation and value (i.e. the "s" and "v"
      in the "HSV" color scheme) to be used for the gradient. The defaults
      are 1.0 for both "sat" and "value".

      In the case of the gray scale gradients, "sat" sets the minimum intensity
      (normally black) and "value" sets the maximum (normally white)

    usage:
      from within PyMOL do "run color_b.py" to load the function definition.
      Then you can use for example:

          color_b (c. a | c. b),mode=ramp,gradient=bwr,nbins=30,sat=.5, value=1.

      to color chains A and B with the Blue-White-Red gradient in 30 colors of equal
      numbers of atoms in each color.
  """

  nbins=int(nbins)
  sat=float(sat)
  value=float(value)
# make sure sat and value are in the range 0-1.0
  sat = min(sat, 1.0)
  sat = max(sat, 0.0)
  value = min(value, 1.0)
  value = max(value, 0.0)
  debug = int(debug)
  if gradient == 'user' and user_rgb == '':
    user_rgb = '50,50,195, 245,245,20, 255,20,20'

# make sure lowercase
  gradient.lower()
  mode.lower()

# Sanity checking
  if nbins == 1:
    #print "\n     WARNING: You specified nbins=1, which doesn't make sense...resetting nbins=11\n"
    nbins=11

  if mode not in ('hist','ramp'):
    #print "\n     WARNING: Unknown mode ",mode, "    ----->   Nothing done.\n"
    return
  elif gradient not in ('bgr','rgb','rainbow','reverserainbow','bwr','rwb','user',
                        'bmr','rmb','rw','wr','ry','yr','gw','wg','bw','wb','gy',
                        'yg','gray','grey','reversegray','reversegrey','user'):
    #print "\n     WARNING: Unknown gradient: ",gradient, "    ----->   Nothing done.\n"
    return

  #print "MODE, GRADIENT, NBINS:", mode,gradient, nbins

# get list of B-factors from selection

# contents of "m.atom[i]": 'b', 'chain', 'coord', 'defaults',
#'elec_radius', 'flags', 'formal_charge', 'get_implicit_valence',
#'get_mass', 'get_number', 'get_signature', 'has', 'hetatm', 'id',
#'in_same_residue', 'index', 'name', 'new_in_residue', 'numeric_type',
#'partial_charge', 'q', 'resi', 'resi_number', 'resn', 'segi', 'ss',
#'stereo', 'symbol', 'u_aniso', 'vdw']
  if cmd.count_atoms(selection) == 0:
    return []
  m = cmd.get_model(selection)
  sel = []
  b_list = []

  if len(m.atom) == 0:
    pass
   # print "Sorry, no atoms selected"

  else:
    if item == 'b':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].b)
    elif item == 'q':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].q)

    elif item == 'partial_charge':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].partial_charge)

    elif item == 'formal_charge':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].formal_charge)

    elif item == 'index':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].index)

    elif item == 'resi':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].resi_number)

    else:
      print("Not configured to work on item %s" % item)
      return

    max_b = float(max(b_list))
    min_b = float(min(b_list))
    #print "Minimum and Maximum B-values: ", min_b, max_b

    if mode == 'ramp':
      # color in bins of equal numbers of atoms
      b_list.sort()

      # subtract 0.1 from the lowest B in order to ensure that the single
      # atom with the lowest B value doesn't get omitted
#      b_list[0] = b_list[0] - 0.1

      bin_num = int(len(b_list)/nbins)
#      sel.append(selection + " and (b < " + str(b_list[bin_num]) + " or b = " + str(b_list[bin_num]) + ")")
      if item == 'index' or item == 'resi':
        sel.append(selection + " and (%s 0-%d" % (item,b_list[bin_num]) + ")")
        for j in range(1,nbins):
          sel.append(selection + " and %s %d-%d" % (item,b_list[j*bin_num],b_list[(j+1)*bin_num]))
      else:
        sel.append(selection + " and (%s < %4.4g" % (item,b_list[bin_num]) + " or %s = %4.4g" % (item,b_list[bin_num]) + ")")
        for j in range(1,nbins):
#        sel.append(selection + " and b > " + str(b_list[j*bin_num]))
          sel.append(selection + " and %s > %4.4g" % (item,b_list[j*bin_num]))
          #print "Color select: ",sel[j]

    elif mode == 'hist':

# check if minimum or maximum was specified and use the entered values
      if minimum != '':
        min_b = float(minimum)
      if maximum != '':
        max_b = float(maximum)
      # histogram:
      # color in bins of equal B-value ranges
      # subtract 0.1 from the lowest B in order to ensure that the single
      # atom with the lowest B value doesn't get omitted
      bin_width = (max_b - min_b)/nbins*1.0
      if item == 'index' or item == 'resi':
        sel.append(selection + " and (%s 0-%d" % (item,min_b + bin_width) + ")")
        for j in range(1,nbins):
          sel.append(selection + " and %s %d-%d" % (item,min_b + j*bin_width,min_b + (j+1)*bin_width))
          #print "Color select: ",sel[j]
      else:
        sel.append(selection + " and (%s < %4.4g" % (item,min_b + bin_width) + " or %s = %4.4g" % (item,min_b + bin_width) + ")")
        for j in range(1,nbins):
          sel.append(selection + " and %s > %4.4g" % (item,min_b + j*bin_width))
          #print "Color select: ",sel[j]

# call the function to create the gradient which returns a list of colours
    colours = make_gradient(sel,gradient,nbins,sat,value, user_rgb,debug)

# do the colouring now
    for j in range(nbins):
      #print "Color select: ",sel[j]
      cmd.color(colours[j],sel[j])
  #    print j,colours[j],sel[j]

def color_q(selection="all",mode="hist",gradient="bgr",nbins=11,sat=1.,value=1.,minimum='',maximum='',debug=0):
  """

  USAGE

    color_q(selection,gradient,mode,nbins,sat,value) ='sel',
      gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb'
      'rw' or 'wr','gw' or 'wg' or 'bw' or 'wb' or 'gy' or 'yg' or 'gray' or 'reversegray' or 'user'
      mode='hist' or 'ramp', q0=0.,q1=1.0,
      nbins=11, sat=1.0, value=1.0)

      This function allows coloring of a selection as a function of
      occupancy.  See color_b for details.
  """
  item='q'
  debug=int(debug)
  color_b(selection,item,mode,gradient,nbins,sat,value,minimum,maximum,debug)

# function for creating the gradient
def make_gradient(sel,gradient,nbins,sat,value,user_rgb,debug=0):
  if gradient == 'bgr' or gradient == 'rainbow':
    gradient = 'bgr'
    col=[]
    coldesc=[]
    for j in range(nbins):
      # must append the str(sel[j]) to the color name so that it is unique
      # for the selection
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # coldesc.append('col' + str(sel[j]) + str(j))

      # create colors using hsv scale (fractional) starting at blue(.6666667)
      # through red(0.00000) in intervals of .6666667/(nbins -1) (the "nbins-1"
      # ensures that the last color is, in fact, red (0)
      # rewrote this to use the colorsys module to convert hsv to rgb
      hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
      #convert to rgb and append to color list
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      #cmd.set_color("col" + gradient + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])
      #print col[j],"defined as ", "col"+str(j)

  elif gradient == 'user':
    # --------------------------------------------
    #  Modified color ramp by Mark Wall 2007.07.20
    # --------------------------------------------
    # assign 3-color ramp values (rgb 0-255)
    # could easily assign RGB color values between 0.0 and 1.0 below
    # !!! Black must be at least 1,0,0 or div by zero error !!!
    #
    #r1, g1, b1 = 255, 255, 225   # low   white
    #r1, g1, b1 = 170, 170, 170   # low   gray
    user_rgb = re.compile('[\[\](){}]').sub('',user_rgb)
    user_rgb_fract = 0
    try:
      r1,g1,b1, r2,g2,b2, r3,g3,b3 = list(map(int,user_rgb.split(',')))
    except ValueError:
      r1,g1,b1, r2,g2,b2, r3,g3,b3 = list(map(float,user_rgb.split(',')))
      user_rgb_fract = 1

    #print 'user_rgb', r1,g1,b1, r2,g2,b2, r3,g3,b3

#    r1, g1, b1 = 50, 50, 195     # low   med blue
#    r2, g2, b2 = 245, 245, 20    # mid   yellow
#    r3, g3, b3 = 255, 20, 20     # high  red
    #
    #r1, g1, b1 = 255, 20, 20     # low   red
    #r2, g2, b2 = 150, 150, 20    # mid   yellow
    #r3, g3, b3 = 20, 20, 195     # high   med blue
    #
    #r1, g1, b1 = 1, 0, 0         # low   black
    #r2, g2, b2 = 155, 155, 155   # mid  gray
    #r3, g3, b3 = 255, 255, 255   # high  white
    #
    #r1, g1, b1 = 0, 50, 200      # low  blue
    #r2, g2, b2 = 1, 0, 0         # mid   black
    #r3, g3, b3 = 255, 255, 20    # high  yellow
    #
    #r1, g1, b1 = 0, 0, 1         # low   black
    #r2, g2, b2 = 200, 0, 0       # mid  red
    #r3, g3, b3 = 255, 255, 0     # high  yellow
    #
    #r1, g1, b1 = 180, 170, 170   # low   gray
    #r2, g2, b2 = 250, 90, 40     # mid  orange
    #r3, g3, b3 = 255, 255, 0     # high  yellow
    #
    #r1, g1, b1 = 235, 255, 255   # low   white
    #r2, g2, b2 = 55, 255, 255    # mid   cyan
    #r3, g3, b3 = 0, 0, 180       # high  blue
    #
    # change color values to fractions
    #
    if max(r1,g1,b1,r2,g2,b2,r3,g3,b3) > 1:
      r1, g1, b1 = float(r1)/255.0, float(g1)/255.0, float(b1)/255.0
      r2, g2, b2 = float(r2)/255.0, float(g2)/255.0, float(b2)/255.0
      r3, g3, b3 = float(r3)/255.0, float(g3)/255.0, float(b3)/255.0

    col=[]
    coldesc=[]
#    print "r1,g1,b1, r2,g2,b2, r3,g3,b3", r1,g1,b1, r2,g2,b2, r3,g3,b3
    for j in range(int(nbins/2)):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from low to mid
      rgb = [r1*((float(nbins)-float(j)*2.0)/float(nbins))+r2*(float(j)*2.0/float(nbins)), \
             g1*((float(nbins)-float(j)*2.0)/float(nbins))+g2*(float(j)*2.0/float(nbins)), \
             b1*((float(nbins)-float(j)*2.0)/float(nbins))+b2*(float(j)*2.0/float(nbins))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
      #print j,"rgb: %4.3f %4.3f %4.3f"% rgb

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

    for j in range(int(nbins/2),nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from mid to high
      rgb = [r2*((float(nbins)-((float(j+1)-float(nbins)/2.0)*2.0))/float(nbins))+r3*(((float(j+1)-float(nbins)/2.0)*2.0)/float(nbins)), \
             g2*((float(nbins)-((float(j+1)-float(nbins)/2.0)*2.0))/float(nbins))+g3*(((float(j+1)-float(nbins)/2.0)*2.0)/float(nbins)), \
             b2*((float(nbins)-((float(j+1)-float(nbins)/2.0)*2.0))/float(nbins))+b3*(((float(j+1)-float(nbins)/2.0)*2.0)/float(nbins))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
      #print j,"rgb: %4.3f %4.3f %4.3f"% rgb

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'rgb' or gradient == 'reverserainbow':
    gradient = 'rgb'
    col=[]
    coldesc=[]
    for j in range(nbins):
      # must append the str(sel[j]) to the color name so that it is unique
      # for the selection
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))

      # create colors using hsv scale (fractional) starting at red(.00000)
      # through blue(0.66667) in intervals of .6666667/(nbins -1) (the "nbins-1"
      # ensures that the last color is, in fact, red (0)
      # rewrote this to use the colorsys module to convert hsv to rgb
      hsv = (colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
      #convert to rgb and append to color list
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'bmr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from blue through magenta to red
      rgb = [min(1.0, float(j)*2/(nbins-1)), 0.0, min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'rmb':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from red through magenta to blue
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), 0.0, min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'rw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from red through white
      rgb = [1.0, float(j)/(nbins-1), float(j)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'wr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from white through red
      rgb = [1.0, float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'ry':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from red through white
      rgb = [1.0, float(j)/(nbins-1), 0]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'yr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from white through red
      rgb = [1.0, float(nbins-j-1)/(nbins-1), 0]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'gw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from green through white
      rgb = [float(j)/(nbins-1), 1.0, float(j)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'wg':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from white through green
      rgb = [float(nbins-j-1)/(nbins-1), 1.0, float(nbins-j-1)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'bw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from blue through white
      rgb = [float(j)/(nbins-1), float(j)/(nbins-1), 1.0 ]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'wb':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from blue through white
      rgb = [float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1), 1.0 ]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'wr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from white through blue
      rgb = [float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1), 1.0]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'gy':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from green through yellow
      rgb = [float(j)/(nbins-1), 1.0, 0.]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'yg':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from green through yellow
      rgb = [float(nbins-j-1)/(nbins-1), 1.0, 0.]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'bwr':
    col=[]
    coldesc=[]
    for j in range(int(nbins/2)):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from blue to white
      rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

    for j in range(int(nbins/2),nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from white to red
      rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'rwb':
    col=[]
    coldesc=[]
    for j in range(int(nbins/2)):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient from red to white
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

    for j in range(int(nbins/2),nbins):
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # coldesc.append('col' + str(sel[j]) + str(j)))
      # create colors in a gradient from white to blue
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'gray' or gradient == 'grey':
# if it is "gray" then sat must be 0!
    sat = 0.0
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient of grays from "sat" to "value"

      hsv = [0, 0, sat + (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

  elif gradient == 'reversegray' or gradient == 'reversegrey':
# if it is "gray" then sat must be 0!
    sat = 0.0
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + gradient + str(j) + str(sel[j]))
      # create colors in a gradient of grays from "sat" to "value"

      hsv = [0, 0, value - (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      if debug:
        print("Colour RGB triplet [ %6.4f, %6.4f, %6.4f ] is defined as %s" % (col[j][0],col[j][1],col[j][2],"col"+str(j)+str(sel[j])))
      cmd.set_color("col" + gradient + str(j) + str(sel[j]),col[j])

#  if debug:
#    for j in range(nbins):
#      print "colour #:",j,"colour RGB triplet: ",col[j]

  #print coldesc
# return the gradient as a list of colors named by their gradient & index (i.e. colbw0,colbw1,colbw2,...)
  return coldesc

# allow calling without parentheses: color_hist_b [selection=], [mode= ],[gradient= ],[nbins= ]
cmd.extend("PDI_Viz_color_b",color_b)
cmd.extend("PDI_Viz_color_q",color_q)
stored.PDI_Viz_color_b = color_b
stored.PDI_Viz_color_q = color_q


