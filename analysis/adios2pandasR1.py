import numpy as np
import adios2
import MDAnalysis as mda
from  MDAnalysis.analysis.rms import RMSD
import sys
import pandas as pd


ref_pdb_file ='../MD_exps/fs-pep/pdb/1FME-0.pdb'

init_pdb = '../MD_exps/fs-pep/pdb/1FME.pdb'

Rs = []
steps = []
fsteps = []
dirs = []

pf = pd.DataFrame(columns=["fstep", "step", "R", "dir"])


with adios2.open("../aggregate/aggregator.bp", "r") as fr:

    nsteps = fr.steps()

    for s in range(nsteps):
    #for s in range(36):
        fstep = next(fr)
        dir = fstep.read('dir').tostring().decode()
        step = fstep.read('step')
        positions = fstep.read('positions')

        outlier_traj = mda.Universe(init_pdb, positions)

        ref_traj = mda.Universe(ref_pdb_file)
        R = RMSD(outlier_traj, ref_traj, select='protein and name CA')
        R.run()

        Rs.append(R.rmsd[:,2][0])
        steps.append(step)
        fsteps.append(s)
        dirs.append(dir)


pf['fstep'] = fsteps
pf['step'] = steps
pf['dir'] = dirs
pf['R'] = Rs

print(pf)

pf.to_csv('pf.csv')



