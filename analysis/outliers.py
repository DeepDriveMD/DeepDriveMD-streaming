import numpy as np
#import adios2
import MDAnalysis as mda
from  MDAnalysis.analysis.rms import RMSD
import sys
import pandas as pd
from multiprocessing import Pool
import glob
import os

ref_pdb_file ='../MD_exps/fs-pep/pdb/1FME-0.pdb'
init_pdb = '../MD_exps/fs-pep/pdb/1FME.pdb'

ref_traj = mda.Universe(ref_pdb_file)

outlier_iterations = map(lambda x: x.split("/")[-1], glob.glob(f"../Outlier_search/archive/*"))

for outlier_iteration in outlier_iterations:
    print(f"outlier_iteration = {outlier_iteration}")
    outlier_fns = glob.glob(f"../Outlier_search/archive/{outlier_iteration}/outlier_pdbs/*.pdb")
    output = []
    for outlier_fn in outlier_fns:
        print(f"outlier_fn = {outlier_fn}")
        outlier_traj = mda.Universe(outlier_fn)
        R = RMSD(outlier_traj, ref_traj, select = 'protein and name CA')
        R.run()
        distance = R.rmsd[:,2][0]
        md5 = os.path.basename(outlier_fn).split(".")[0]
        output.append((md5, distance))

    df = pd.DataFrame(output, columns=["md5","RMSD"])

    print(f"number of outliers={len(output)}")

    df.to_csv(f"outliers_{outlier_iteration}.csv")

