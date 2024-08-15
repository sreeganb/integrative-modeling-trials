#!/usr/bin/env python

import numpy as np

# Load PDB files 
load "glut9.pdb", reference 
load "glut9.B99990001.pdb", object2 
#load protein3.pdb, object3 

# Align objects object2 and object3 to object1 
align object2, reference 
#align object3, reference

