# Modeller script for QTY mutations
from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients

log.verbose()
env = Environ()
env.io.atom_files_directory = ['.']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read the original model
mdl = Model(env, file='/Users/aneeshbondugula/Downloads/QTY_docking/01_prepare_structure/output/8EF5_receptor_ligand_fixed.pdb')

# Define mutations
mutations = [(6, 'T'), (8, 'T'), (11, 'Q'), (14, 'T'), (15, 'T'), (17, 'T'), (18, 'T'), (20, 'Q'), (21, 'Y'), (24, 'Y'), (25, 'Q'), (26, 'T'), (42, 'T'), (44, 'T'), (45, 'Y'), (47, 'Q'), (49, 'Q'), (53, 'Q'), (58, 'Q'), (60, 'Y'), (63, 'T'), (66, 'Q'), (79, 'T'), (80, 'T'), (81, 'T'), (83, 'T'), (89, 'Y'), (92, 'T'), (93, 'Y'), (95, 'Q'), (100, 'T'), (123, 'T'), (124, 'T'), (126, 'T'), (130, 'T'), (131, 'Q'), (135, 'T'), (137, 'Q'), (139, 'T'), (141, 'Y'), (168, 'Q'), (169, 'Q'), (171, 'T'), (173, 'T'), (174, 'Y'), (175, 'T'), (176, 'Y'), (178, 'Y'), (179, 'T'), (182, 'T'), (183, 'Q'), (184, 'T'), (185, 'T'), (187, 'T'), (215, 'T'), (219, 'T'), (220, 'Q'), (221, 'T'), (222, 'T'), (223, 'T'), (225, 'T'), (226, 'Y'), (227, 'T'), (228, 'T'), (233, 'T'), (235, 'T'), (237, 'T'), (238, 'T'), (239, 'T'), (250, 'Y'), (253, 'T'), (257, 'Y'), (259, 'T'), (261, 'Q'), (268, 'Q'), (271, 'T'), (272, 'Q')]

# Apply mutations
for pos, new_res in mutations:
    sel = Selection(mdl.residues[str(pos) + ':' + 'A'])
    sel.mutate(residue_type=new_res)

# Optimize the mutated sidechains
mdl.write(file='02_introduce_mutations/output/qty_mutant.pdb')
