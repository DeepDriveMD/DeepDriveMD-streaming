import os, random, json, shutil, time, subprocess, sys
import argparse 
import numpy as np 
from glob import glob
import MDAnalysis as mda
from utils import outliers_from_cvae, adios_to_cvae, adios_list_to_cvae  
from utils import predict_from_cvae, outliers_from_latent
from utils import write_pdb_frame, make_dir_p 
from  MDAnalysis.analysis.rms import RMSD
import hashlib
import pickle
#from OutlierDB import *
#import mytimer
from lockfile import LockFile

lastN = int(sys.argv[1])

model_best = os.path.basename(os.path.realpath('best_model'))

model_dim = int(model_best[10:12]) 
model_best = f'{model_best}/best.h5'
adios_files_list = glob("../aggregate/aggregator*.bp")

cvae_input = adios_list_to_cvae(adios_files_list, lastN)

print(cvae_input[0].shape)

cm_predict = predict_from_cvae(model_best, cvae_input[0], hyper_dim=model_dim) 

np.save("projection.npy", cm_predict)
