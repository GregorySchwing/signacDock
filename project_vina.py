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
import random
from vina import Vina
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
        trajFile.reset()
        trajFile.setAtoms(mdm2.protein)
        writeDCD(join(jobprefix,'mdm2_merged_protein.dcd'), trajFile)
        writePDB(join(jobprefix,'ref_protein.pdb'), mdm2.protein)
for ph in phs:
    jobprefix = os.path.join(dcd_prefix, str(ph))
    refjobprefix = os.path.join(dcd_prefix, str(5.63))
    if not(os.path.isfile(os.path.join(jobprefix,"mdm2_bb_aligned.dcd"))):
        pdbFile = join(refjobprefix, 'ref.pdb')
        print(pdbFile)
        mdm2_bb = parsePDB(pdbFile)
        trajFile = Trajectory(join(jobprefix,'mdm2_merged_bb.dcd'))
        trajFile.link(mdm2_bb)
        out_bb = DCDFile(join(jobprefix,'mdm2_bb_aligned.dcd'), 'w')
        for frame in trajFile:
            frame.superpose()
            out_bb.write(mdm2_bb)
    if not(os.path.isfile(os.path.join(jobprefix,"mdm2_protein_aligned.dcd"))):
        pdbFile = join(jobprefix, 'ref_protein.pdb')
        print(pdbFile)
        mdm2 = parsePDB(pdbFile)
        trajFile = Trajectory(join(jobprefix,'mdm2_merged_protein.dcd'))
        trajFile.link(mdm2)
        out_bb = DCDFile(join(jobprefix,'mdm2_merged_protein_aligned.dcd'), 'w')
        for frame in trajFile:
            frame.superpose()
            out_bb.write(mdm2)
vina = False    

for ph in phs:
    jobprefix = os.path.join(dcd_prefix, str(ph))
    if (os.path.isfile(os.path.join(jobprefix,"mdm2_merged_protein_aligned.dcd"))):
        jobprefix = os.path.join(dcd_prefix, str(ph))
        pdbFile = join(jobprefix, 'ref_protein.pdb')
        mdm2_bb = parsePDB(pdbFile)
        trajFile = Trajectory(join(jobprefix,'mdm2_merged_protein_aligned.dcd'))
        print(pdbFile)
        trajFile.link(mdm2_bb)
        random.seed(10)
        print(random.random()) 
        for N in [1, 10, 100, 500]:
            # Generate 10 unique random numbers within a range
            print(trajFile.numFrames())
            num_list = random.sample(range(0, trajFile.numFrames()), N)
            for n in range(1, N+1, 1):
                loc = "view/pH/%s/N/%s/n/%s" % (ph, N, n)
                frame = join(dcd_prefix,loc)
                writePDB(join(frame,"receptor.pdb"), trajFile[num_list[n-1]])
                if (vina):
                    inputFile = join(frame,"receptor.pdb")
                    outputFile = join(frame,"receptor.pdbqt")
                    cmd = "prepare_receptor -r {inf} -o {of}".format(inf = inputFile, of = outputFile)
                    subprocess.run(cmd, shell=True, check=True)
for N in [1, 10, 100]:
    for ph in phs:
        jobprefix = os.path.join(dcd_prefix, str(ph))
        ligands = os.path.join(dcd_prefix, "PossessActVal.sdf")
        # Generate 10 unique random numbers within a range
        trajFile = Trajectory(join(jobprefix,'mdm2_merged_protein_aligned.dcd'))
        print(trajFile.numFrames())
        num_list = random.sample(range(0, trajFile.numFrames()), N)
        for n in range(1, N+1, 1):            
            loc = "view/pH/%s/N/%s/n/%s" % (ph, N, n)
            frame = join(dcd_prefix,loc)
            recFilePDB = join(frame,"receptor.pdb")
            mdm2 = parsePDB(recFilePDB)
            bsatoms = mdm2.select('resnum 913:925')
            com = calcCenter(bsatoms)
            recFile = join(frame,"receptor.pdbqt")
            #dockFile = join(frame,"bs_HL1.sdf.gz")
            #dockFile = join(frame,"whole_docked.sdf.gz")
            # whole protein
            #cmd = "gnina -r {recf} -l {ligf} --autobox_ligand {recf} -o {docf} --exhaustiveness 64 --verbosity 4".format(recf = recFile, ligf = ligands, docf = dockFile)
            v.dock(exhaustiveness=32, n_poses=20)
            vina --receptor 1iep_receptor.pdbqt --ligand 1iep_ligand.pdbqt \
       --config 1iep_receptor_vina_box.txt \
       --exhaustiveness=32 --out 1iep_ligand_vina_out.pdbqt
            cmd = "gnina -r {recf} -l {ligf} --center_x {cenx} --center_y {ceny} --center_z {cenz} --size_x 20 --size_y 20 --size_z 20 -o {docf} --exhaustiveness 8 --verbosity 4".format(recf = recFile, ligf = ligands, docf = dockFile, cenx = com[0], ceny = com[1], cenz = com[2],)
            print(cmd)
            quit()
            #subprocess.run(cmd, shell=True, check=True)      

@FlowProject.label
def volume_computed(job):
    return job.isfile("volume.txt")
    

@volume_group
@FlowProject.operation
@FlowProject.post(volume_computed)
def compute_volume(job):
    volume = job.sp.N
    with open(job.fn('volume.txt'), 'w') as file:
        file.write(str(volume) + '\n')
        
@volume_group
@FlowProject.operation
@FlowProject.pre(volume_computed)
@FlowProject.post.isfile("data.json")
def store_volume_in_json_file(job):
    with open(job.fn("volume.txt")) as textfile:
        data = {"volume": float(textfile.read())}

        with open(job.fn("data.json"), "w") as jsonfile:
            json.dump(data, jsonfile)


@volume_group
@FlowProject.operation
@FlowProject.pre.after(compute_volume)
@FlowProject.post(lambda job: 'volume' in job.document)
def store_volume_in_document(job):
    with open(job.fn("volume.txt")) as textfile:
        job.document.volume = float(textfile.read())
        




if __name__ == "__main__":
    FlowProject().main()
    
    
