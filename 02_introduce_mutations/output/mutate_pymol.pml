# PyMOL Mutation Script for QTY μOR
# Run with: pymol -cq mutation_script.pml

# Load structure
load /Users/aneeshbondugula/Downloads/QTY_docking/01_prepare_structure/output/8EF5_receptor_ligand_fixed.pdb, receptor

# Mutations to apply

# Mutation: I6T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//6/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I8T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//8/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L11Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//11/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I14T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//14/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V15T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//15/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V17T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//17/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V18T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//18/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L20Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//20/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F21Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//21/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F24Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//24/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L25Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//25/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V26T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//26/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I42T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//42/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I44T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//44/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F45Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//45/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L47Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//47/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L49Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//49/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L53Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//53/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L58Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//58/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F60Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//60/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V63T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//63/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L66Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//66/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I79T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//79/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V80T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//80/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I81T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//81/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I83T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//83/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F89Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//89/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I92T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//92/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F93Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//93/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L95Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//95/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V100T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//100/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I123T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//123/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I124T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//124/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V126T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//126/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I130T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//130/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L131Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//131/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I135T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//135/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L137Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//137/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V139T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//139/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F141Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//141/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L168Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//168/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L169Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//169/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I171T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//171/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V173T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//173/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F174Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//174/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I175T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//175/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F176Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//176/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F178Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//178/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I179T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//179/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V182T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//182/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L183Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//183/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I184T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//184/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I185T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//185/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V187T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//187/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I215T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//215/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V219T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//219/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L220Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//220/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V221T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//221/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V222T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//222/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V223T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//223/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V225T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//225/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F226Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//226/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I227T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//227/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V228T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//228/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I233T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//233/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I235T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//235/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V237T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//237/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I238T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//238/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I239T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//239/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F250Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//250/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V253T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//253/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: F257Y
wizard mutagenesis
cmd.get_wizard().set_mode("TYR")
cmd.get_wizard().do_select("/A//257/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: I259T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//259/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L261Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//261/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L268Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//268/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: V271T
wizard mutagenesis
cmd.get_wizard().set_mode("THR")
cmd.get_wizard().do_select("/A//271/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Mutation: L272Q
wizard mutagenesis
cmd.get_wizard().set_mode("GLN")
cmd.get_wizard().do_select("/A//272/")
cmd.get_wizard().apply()
cmd.set_wizard()

# Save mutated structure
save 02_introduce_mutations/output/qty_mutant.pdb, receptor

# Quit
quit
