import os, random, json, shutil, time, subprocess, sys
import numpy as np 
from glob import glob
import MDAnalysis as mda
from utils import outliers_from_cvae, adios_to_cvae, adios_list_to_cvae  
from  MDAnalysis.analysis.rms import RMSD
from multiprocessing import Pool

ps = int(sys.argv[1])
lastN = int(sys.argv[2])

ref_pdb = '../MD_exps/fs-pep/pdb/1FME-0.pdb'
init_pdb = '../MD_exps/fs-pep/pdb/1FME.pdb'

def f(position):
    outlier_traj = mda.Universe(init_pdb, position)
    ref_traj = mda.Universe(ref_pdb)
    R = RMSD(outlier_traj, ref_traj, select = 'protein and name CA')
    R.run()
    return R.rmsd[:,2][0]

model_best = os.path.basename(os.path.realpath('best_model'))

model_dim = int(model_best[10:12]) 
model_best = f'{model_best}/best.h5'
adios_files_list = glob("../aggregate/aggregator*.bp")

cvae_input = adios_list_to_cvae(adios_files_list, lastN)

positions = cvae_input[1]

with Pool(processes=ps) as pool:
    Rs = pool.map(f, positions)

np.save('rmsd.npy', np.array(Rs))
