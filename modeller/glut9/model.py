#!/usr/bin/env/python

# Comparative modeling by the AutoModel class
from modeller import * # Load standard Modeller classes
from modeller.automodel import * # Load the AutoModel class
log.verbose() # request verbose output
env = Environ() # create a new MODELLER environment to build this model in
#-------------------------------------------------
# Create the alignment files for modeller to use
#-------------------------------------------------
aln = Alignment(env)
mdl = Model(env, file='glut3', model_segment=('FIRST:A', 'LAST:A'))
aln.append_model(mdl, align_codes='glut3', atom_files='glut3.pdb')
aln.append(file='glut9.fasta', align_codes='glut9', alignment_format='FASTA')
aln.align2d(max_gap_length=50)
aln.write(file='Target_template.ali', alignment_format='PIR')


## directories for input atom files
#env.io.atom_files_directory = [’.’, ’../atom_files’]
#a = AutoModel(env,
#alnfile = ’alignment.ali’, # alignment filename
#knowns = ’5fd1’, # codes of the templates
#sequence = ’1fdx’) # code of the target
#a.starting_model= 1 # index of the first model
#a.ending_model = 1 # index of the last model
## (determines how many models to calculate)
#a.make() # do the actual comparative modeling
