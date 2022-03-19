# project.py
from flow import FlowProject
import json
import signac
import os
import glob
from os.path import join
from natsort import natsorted # pip install natsort
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import AnalysisFromFunction
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader
import subprocess
from prody import *
import matplotlib.pyplot as mpl

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

volume_group = FlowProject.make_group(name='volume')
dcd_prefix = signac.get_project().root_directory()

# Move this into a post condition of a namd run job
phs = [5.63, 7.1]

for ph in phs:
    jobprefix = os.path.join(dcd_prefix, str(ph))
    if not(os.path.isfile(os.path.join(jobprefix,"mdm2_merged_bb.dcd"))):
        pdbFile = join(jobprefix, 'step3_input.pdb')
        print(pdbFile)
        mdm2 = parsePDB(pdbFile)
        trajfiles = glob.glob(join(jobprefix, '*.dcd'))
        trajfilessorted = natsorted(trajfiles)
        trajFile = Trajectory(trajfilessorted[0])
        for traj in trajfilessorted[1:]:
            trajFile.addFile(traj)
        trajFile.link(mdm2)
        trajFile.setAtoms(mdm2.backbone)
        writeDCD(join(jobprefix,'mdm2_merged_bb.dcd'), trajFile)
        writePDB(join(jobprefix,'ref.pdb'), mdm2.backbone)
        trajFile.setAtoms(mdm2.noh)
        writeDCD(join(jobprefix,'mdm2_merged_bb_noh.dcd'), trajFile)
        writePDB(join(jobprefix,'ref_noh.pdb'), mdm2.noh)
for ph in phs:
    jobprefix = os.path.join(dcd_prefix, str(ph))
    refjobprefix = os.path.join(dcd_prefix, str(5.63))
    if not(os.path.isfile(os.path.join(jobprefix,"mdm2_merged_bb_aligned.dcd"))):
        pdbFile = join(refjobprefix, 'ref.pdb')
        print(pdbFile)
        mdm2_bb = parsePDB(pdbFile)
        trajFile = Trajectory(join(jobprefix,'mdm2_merged_bb.dcd'))
        trajFile.link(mdm2_bb)
        out_bb = DCDFile(join(jobprefix,'mdm2_bb_aligned.dcd'), 'w')
        for frame in trajFile:
            frame.superpose()
            out_bb.write(mdm2_bb)
    if not(os.path.isfile(os.path.join(jobprefix,"mdm2_merged_protein_aligned.dcd"))):
        pdbFile = join(refjobprefix, 'ref_noh.pdb')
        print(pdbFile)
        mdm2_bb = parsePDB(pdbFile)
        trajFile = Trajectory(join(jobprefix,'mdm2_merged_protein.dcd'))
        trajFile.link(mdm2_bb)
        out_bb = DCDFile(join(jobprefix,'mdm2_merged_protein_aligned.dcd'), 'w')
        for frame in trajFile:
            frame.superpose()
            out_bb.write(mdm2_bb)
            
refjobprefix = os.path.join(dcd_prefix, str(5.63))            
structure = parsePDB(join(refjobprefix, 'ref.pdb'))
trajectory = Trajectory(join(refjobprefix,"mdm2_bb_aligned.dcd"))
jobprefix = os.path.join(dcd_prefix, str(7.1))
trajectory.addFile(join(jobprefix,'mdm2_bb_aligned.dcd'))

trajectory.link(structure)
trajectory.setCoords(structure)
trajectory.setAtoms(structure.calpha)
eda = EDA('mdm2')
eda.buildCovariance( trajectory )
eda.calcModes()
saveModel(eda)

for mode in eda[:4]:
    print(calcFractVariance(mode).round(2))


mdm2ca_sim1 = trajectory[:500]
mdm2ca_sim1.superpose()
mdm2ca_sim2 = trajectory[500:]
mdm2ca_sim2.superpose()

showProjection(mdm2ca_sim1, eda[:3], color='red', marker='.');
showProjection(mdm2ca_sim2, eda[:3], color='blue', marker='.');

writeNMD('mdm2_eda.nmd', eda[:3], structure.select('calpha'))
mpl.show()
import csv
import numpy as np
refjobprefix = os.path.join(dcd_prefix, str(5.63))  
refpdbFile = join(refjobprefix, 'ref.pdb')          
mdm2 = parsePDB(pdbFile)
resnums = mdm2.getResnums()
resnumsUn = f7(resnums)
arrays = []  
for ph in phs:
    jobprefix = os.path.join(dcd_prefix, str(ph))
    traj = join(jobprefix,'mdm2_bb_aligned.dcd')
    ensemble = parseDCD(traj)
    ensemble.setAtoms(mdm2.backbone)
    ensemble.setCoords(mdm2.backbone)
    ensemble.superpose()
    rmsf = ensemble.getRMSFs()
    arrays.append(rmsf)
    wtr = csv.writer(open (join(jobprefix,"rmsf.csv"), 'w'), delimiter=',', lineterminator='\n')
    for x, y in zip(resnumsUn, rmsf) : wtr.writerow ([x, y])
    
net = np.subtract(arrays[0], arrays[1])
wtr = csv.writer(open ("netrmsf.csv", 'w'), delimiter=',', lineterminator='\n')
for x, y in zip(resnumsUn, net): wtr.writerow ([x, y])

wtr = csv.writer(open ("netrmsfjustrmsf.csv", 'w'), delimiter=',', lineterminator='\n')
for x, y in zip(resnumsUn, net): wtr.writerow ([y])


from MDAnalysis import Universe

import MDAnalysis.analysis.encore as encore

#from MDAnalysis.tests.datafiles import PDB, DCD, DCD2


ens1 = Universe("5.63/ref.pdb", "5.63/mdm2_bb_aligned.dcd")

ens2 = Universe("5.63/ref.pdb", "7.1/mdm2_bb_aligned.dcd")
# ARG 366 = ARG 913
# ARG 367 = ARG 914
# LYS 371 = LYS 918
# ARG 372 = ARG 919
# PHE 374 = PHE 921
# PHE 375 = PHE 922
# LEU 377 = LEU 924
# ARG 378 = ARG 925
HES, details = encore.hes([ens1, ens2], select="resid 913:914 or resid  918 or resid 919:922 or resid 924:925")
print(HES)
print(details)
# whole protein
#DRES, details = encore.dres([ens1,ens2])
# JKY binding site
DRES, details = encore.dres([ens1,ens2], estimate_error=True, bootstrapping_samples=20, calc_diagonal=True, ncores=10, select="resid 913:914 or resid  918 or resid 919:922 or resid 924:925")
print("DRES")
print(DRES)
print("error")
print(details)
wtr = csv.writer(open ("6.error.csv", 'w'), delimiter=',', lineterminator='\n')
for x in details: wtr.writerow ([x])
# whole protein
#DRES, details = encore.dres([ens1,ens2])
# JKY binding site
DRES, details = encore.dres([ens1,ens2], select="resid 913:914 or resid  918 or resid 919:922 or resid 924:925")
print("DRES")
print(DRES)
print("error")
print(details)
redCoords = details['reduced_coordinates']
conMem = details['dimensionality_reduction_details']

print(redCoords)
print(len(redCoords))
print(conMem)
print(len(conMem))
memArray = conMem['ensemble_membership']
print(memArray)
print(len(memArray))

fig = mpl.figure()
ax = fig.add_subplot(projection='3d')
scatter = ax.scatter(redCoords[0][0], redCoords[0][1], redCoords[0][2],
                color=[["red", "blue"][m-1] for m
                in memArray])
ax.set_title("Myc-Max Titration State Ensembles", pad=1.0)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
mpl.show()

dres_conv1 = encore.dres_convergence(ens1, 10)
wtr = csv.writer(open ("6.53_convergence.csv", 'w'), delimiter=',', lineterminator='\n')
for x in dres_conv1: wtr.writerow (x)
dres_conv2 = encore.dres_convergence(ens2, 10)
wtr = csv.writer(open ("7.1_convergence.csv", 'w'), delimiter=',', lineterminator='\n')
for x in dres_conv2: wtr.writerow (x)

