import numpy as np
import adios2
import MDAnalysis as mda
from  MDAnalysis.analysis.rms import RMSD
import sys
import pandas as pd
from multiprocessing import Pool
import glob

ref_pdb_file ='../MD_exps/fs-pep/pdb/1FME-0.pdb'

init_pdb = '../MD_exps/fs-pep/pdb/1FME.pdb'

Rs = []
steps = []
fsteps = []
dirs = []

pf = pd.DataFrame(columns=["fstep", "step", "R", "dir"])


bpfiles = glob.glob("../aggregate/*.bp")


STEPS = []
POSITIONS = []

for bp in bpfiles:
    with adios2.open(bp, "r") as fr:
        n = fr.steps()
        vars = fr.available_variables()
        
        print(vars)

        name = 'positions'
        shape = list(map(int, vars[name]['Shape'].split(",")))
        zs = list(np.zeros(len(shape), dtype=np.int))
        positions = fr.read(name, zs, shape, 0, n)
        print(type(positions))
        print(positions.shape)
        sys.stdout.flush()
        POSITIONS.append(positions)
        
        name = 'step'
        steps = fr.read(name, [],[], 0, n)
        print(type(steps))
        print(steps.shape)
        print(steps)
        STEPS.append(steps)

positions = np.concatenate(POSITIONS)
steps = np.concatenate(STEPS)

def f(position):
    outlier_traj = mda.Universe(init_pdb, position)
    ref_traj = mda.Universe(ref_pdb_file)
    R = RMSD(outlier_traj, ref_traj, select = 'protein and name CA')
    R.run()
    return R.rmsd[:,2][0]


with Pool(processes=42) as pool:
    Rs = pool.map(f, positions)

print(type(Rs))

pf['fstep'] = list(np.arange(len(steps)))
pf['step'] = steps
#pf['dir'] = dirs
pf['R'] = Rs

print(pf)

pf.to_csv('pf.csv')



'''




with adios2.open("../aggregate/aggregator.bp", "r") as fr:
    positions = fr.read('positions')
    print(len(positions))
    print(positions.shape)

sys.exit(0)


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

'''

