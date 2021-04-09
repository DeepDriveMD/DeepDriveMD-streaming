import simtk.unit as u
import sys, os, shutil 
import argparse 
import time
import subprocess
from OutlierDB import *
import pickle
import mytimer
from lockfile import LockFile

from MD_utils.openmm_simulation import openmm_simulate_amber_fs_pep 

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--pdb_file", dest="f", help="pdb file")
parser.add_argument("-p", "--topol", dest='p', help="topology file")
parser.add_argument("-c", help="check point file to restart simulation")
parser.add_argument("-l", "--length", default=50, help="how long (ns) the system will be simulated")
parser.add_argument("-g", "--gpu", default=0, help="id of gpu to use for the simulation")
parser.add_argument("-d", "--dir", dest="d", help="top directory")
parser.add_argument("-i", "--init", dest="i", action="store_true")
parser.add_argument("-t", "--task", dest="t", help="taskid to use in naming directories to avoid race condition")
args = parser.parse_args() 

if args.f: 
    pdb_file = os.path.abspath(args.f) 
else: 
    raise IOError("No pdb file assigned...") 

if args.p: 
    top_file = os.path.abspath(args.p) 
else: 
    top_file = None 

if args.c: 
    check_point = os.path.abspath(args.c) 
else: 
    check_point = None 

if args.d:
    top_dir = os.path.abspath(args.d)
else:
    top_dir = None

if args.t:
    task_id = args.t
else:
    task_id = None


# pdb_file = os.path.abspath('./pdb/100-fs-peptide-400K.pdb')
# ref_pdb_file = os.path.abspath('./pdb/fs-peptide.pdb')

gpu_index = 0 # os.environ["CUDA_VISIBLE_DEVICES"]

# check_point = None


def find_pdb(top_dir):
    dbfn = f'{top_dir}/Outlier_search/OutlierDB.pickle'
    mylock1 = LockFile(dbfn)
    if(os.path.exists(dbfn)):
        # my_lock(dbfn)
        mytimer.mytime_label("lock_wait",start=1)
        mylock1.acquire()
        mytimer.mytime_label("lock_wait",start=-1)
        with open(dbfn, 'rb') as f:
            db = pickle.load(f)
        pdb_file = db.next()

        with open(dbfn, 'wb') as f:
            pickle.dump(db, f)

        # my_unlock(dbfn)
        mylock1.release()
    else:
        pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/1FME.pdb'
        #pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/best_12h.pdb'
    return pdb_file


def find_random_pdb(top_dir):
    dbfn = f'{top_dir}/Outlier_search/OutlierDB.pickle'
    mylock1 = LockFile(dbfn)
    if(os.path.exists(dbfn)):
        mytimer.mytime_label("lock_wait",start=1)
        mylock1.acquire()
        mytimer.mytime_label("lock_wait",start=-1)
        with open(dbfn, 'rb') as f:
            db = pickle.load(f)
        mylock1.release()
        pdb_file = db.next_random()
    else:
        #pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/100-fs-peptide-400K.pdb'
        pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/1FME.pdb'
        #pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/best_12h.pdb'
    return pdb_file


def find_best_pdb(top_dir):
    dbfn = f'{top_dir}/Outlier_search/OutlierDB.pickle'
    mylock1 = LockFile(dbfn)
    if(os.path.exists(dbfn)):
        mytimer.mytime_label("lock_wait",start=1)
        mylock1.acquire()
        mytimer.mytime_label("lock_wait",start=-1)
        with open(dbfn, 'rb') as f:
            db = pickle.load(f)
        mylock1.release()
        pdb_file = db.next_10best()
    else:
        #pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/100-fs-peptide-400K.pdb'
        pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/1FME.pdb'
        #pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/best_12h.pdb'
    return pdb_file
    
def prepare(top_dir):
    print(f"in prepare top_dir={top_dir}")
    os.chdir(f'{top_dir}/MD_exps/fs-pep')
    print(subprocess.getstatusoutput("mkdir -p new all running stopped"))
    os.chdir('all')

    time_stamp = time.time()
    run_dir = f'omm_runs_{time_stamp}_{task_id}'
    mytimer.mytime_label("lock4", start=1)
    while(os.path.exists(run_dir)):
        time.sleep(0.1)
        time_stamp = time.time()
        run_dir = f'omm_runs_{time_stamp}_{task_id}'
    mytimer.mytime_label("lock4", start=-1)
    os.mkdir(run_dir)
    os.chdir(run_dir)
    subprocess.getstatusoutput(f'cp {top_dir}/adios.xml .')
    f = open('adios.xml','r')
    textxml = f.read()
    f.close()
    textxml = textxml.replace("SimulationOutput",run_dir)
    print(textxml)

    f = open('adios.xml', 'w')
    f.write(textxml)
    f.close()

    dir_in_all = os.getcwd()
    link_in_new = dir_in_all.replace("all","new")
    subprocess.getstatusoutput(f"ln -s {dir_in_all} {link_in_new}")
    if(not args.i):
        subprocess.getstatusoutput(f"cp {pdb_file} ./")
    sys.stdout.flush()
    sys.stderr.flush()
    return dir_in_all


print("hostname = ", subprocess.getstatusoutput("hostname")[1])

while(not os.path.exists(f'{top_dir}/aggregate/stop.aggregator')):
    mytimer.mytime_label("prepare", start=1)
    run_dir = prepare(top_dir)
    mytimer.mytime_label("prepare", start=-1)
    mytimer.mytime_label("find_best_pdb", start=1)
    pdb_file = find_best_pdb(top_dir)
    mytimer.mytime_label("find_best_pdb", start=-1)
    print(f"As initial condition using {pdb_file}")
    subprocess.getstatusoutput(f"touch {run_dir}/start.simulation")
    print(f"length = {args.length}")
    print(f"Running in {run_dir}"); sys.stdout.flush()
    mytimer.mytime_label("openmm_simulate", start=1)
    openmm_simulate_amber_fs_pep(pdb_file,
                                 check_point = check_point,
                                 GPU_index=gpu_index,
                                 output_traj="output.dcd",
                                 output_log="output.log",
                                 output_cm='output_cm.h5',
                                 report_time=50*u.picoseconds,
                                 sim_time=float(args.length)*u.nanoseconds)
    mytimer.mytime_label("openmm_simulate", start=-1)    

