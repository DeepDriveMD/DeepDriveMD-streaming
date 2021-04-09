#from threadpoolctl import threadpool_info, threadpool_limits
import os, random, json, shutil, time, subprocess, sys
import argparse 
from pprint import pprint
import numpy as np 
from glob import glob
import MDAnalysis as mda
from utils import outliers_from_cvae, adios_to_cvae, adios_list_to_cvae  
from utils import predict_from_cvae, outliers_from_latent
from utils import write_pdb_frame, make_dir_p 
from  MDAnalysis.analysis.rms import RMSD
import hashlib
import pickle
from OutlierDB import *
import mytimer
from lockfile import LockFile
from aggregator_reader import *

from dask.distributed import Client, wait

def archive_outliers(j):
    archive_dir = f"../archive/{j}"
    subprocess.getstatusoutput(f"mkdir -p {archive_dir}")
    print(subprocess.getstatusoutput(f"cp *.json {archive_dir}"))
    print(subprocess.getstatusoutput(f"mkdir -p {archive_dir}/outlier_pdbs"))
    print(subprocess.getstatusoutput(f"cp outlier_pdbs/*.pdb {archive_dir}/outlier_pdbs/"))
    print(subprocess.getstatusoutput(f"cp *.pickle {archive_dir}/"))

def write_pdb(myframe, hash, myframe_v, pdb_file, outliers_pdb_path):
    outlier_pdb_file = f'{outliers_pdb_path}/{hash}.pdb'
    outlier_v_file = f'{outliers_pdb_path}/{hash}.npy'
    write_pdb_frame(myframe, pdb_file, outlier_pdb_file)
    np.save(outlier_v_file, myframe_v)
    return 0

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--md", help="Aggreggated directory")
    parser.add_argument("-c", "--cvae", help="Input: CVAE model directory")
    parser.add_argument("-p", "--pdb", help="Input: pdb file") 
    parser.add_argument("-r", "--ref", default=None, help="Input: Reference pdb for RMSD") 
    parser.add_argument("-n", "--nbp_files", dest="n", default=1, type=int, help="Number of aggregated bpfiles")
    args = parser.parse_args()
    return args

def wait_for_bp(md, n, waiting_time1):
    while(True):
        adios_files_list = mytimer.get_bp_files(md, n)
        if(adios_files_list == None):
            print(f"Waiting for the aggregated bp file to be produced")
            sys.stdout.flush()
            time.sleep(waiting_time1)
        else:
            break
    time.sleep(waiting_time1)
    return adios_files_list

def wait_for_bp_nonempty(adios_files_list, waiting_time1):
    while(True):
        try:
            steps = mytimer.bp_steps(adios_files_list)
            break
        except Exception as e:
            print(e)
            print("Waiting for the aggregated bp file to have data")
            time.sleep(waiting_time1)
    return steps

def wait_for_enough_data(adios_files_list, steps, waiting_time2, current_step, min_step_increment):
    while(steps - current_step < min_step_increment):
        print(f"Waiting for enough time steps to accumulate since the last time: {steps} - {current_step} < {min_step_increment}")
        sys.stdout.flush()
        time.sleep(waiting_time2)
        steps = mytimer.bp_steps(adios_files_list)
    return steps

def wait_for_model(cvae, waiting_time3):
    model_weights = sorted(glob(os.path.join(args.cvae, 'cvae_runs_*/best.h5')))
    while(len(model_weights) == 0):
        print(f"Waiting for the model to become available")
        sys.stdout.flush()
        time.sleep(waiting_time3)
        model_weights = sorted(glob(os.path.join(cvae, 'cvae_runs_*/best.h5')))

    model_best = model_weights[0]
    loss_model_best = np.load(os.path.join(os.path.dirname(model_best), 'loss.npy'))[-1]
    print(f'type(loss_model_best) = {type(loss_model_best)}')
    print(loss_model_best)

    for i in range(len(model_weights)):
        if i + 1 < len(model_weights):
            if int(os.path.basename(os.path.dirname(model_weights[i]))[10:12]) != \
               int(os.path.basename(os.path.dirname(model_weights[i+1]))[10:12]):
                loss = np.load(os.path.join(os.path.dirname(model_weights[i]), 'loss.npy'))[-1]
                if loss < loss_model_best:
                    model_best, loss_model_best = model_weights[i], loss
        else:
            loss = np.load(os.path.join(os.path.dirname(model_weights[i]), 'loss.npy'))[-1]
            if loss < loss_model_best:
                model_best, loss_model_best = model_weights[i], loss
        print(f"i = {i}, model_best = {model_best}, loss_model_best = {loss_model_best}")
    print( "Using model {} with loss {}".format(model_best, loss_model_best) )
    sys.stdout.flush()

    return model_best

def run_dir():
    top_dir = os.getcwd()
    print(f"top_dir = {top_dir}")
    tmp_dir = 'tmp'
    if(not os.path.exists(tmp_dir)):
        os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    return top_dir, tmp_dir

def outlier_read(mystreams):
    mytimer.mytime_label("outlier_read", 1)
    cvae_input = mystreams.next()
    mytimer.mytime_label("outlier_read", -1)
    for ci in cvae_input:
        print(f"shape = {ci.shape}")
    return cvae_input

def eps_init(model_best):
    eps_record_filepath = '../eps_record.json' 
    if os.path.exists(eps_record_filepath): 
        eps_file = open(eps_record_filepath, 'r')
        eps_record = json.load(eps_file) 
        eps_file.close() 
    else: 
        eps_record = {} 

    model_dim = int(os.path.basename(os.path.dirname(model_best))[10:12]) 
    print( 'Model latent dimension: %d' % model_dim  )

    if str(model_best) in eps_record.keys(): 
        eps = eps_record[model_best] 
    else: 
        eps = 1.3

    return eps, eps_record, model_dim, eps_record_filepath

def predict(model_best, cvae_input, model_dim):
    mytimer.mytime_label("outlier_predict", 1)
    cm_predict = predict_from_cvae(model_best, cvae_input[0], hyper_dim=model_dim) 
    mytimer.mytime_label("outlier_predict", -1)
    return cm_predict

def cluster(cm_predict, outlier_list, eps_record, eps, model_best,
            outlier_count, outlier_max, outlier_min):
    mytimer.mytime_label("outlier_eps", 1)
    
    while outlier_count > 0:
        n_outlier = 0
        try:
            outliers = np.squeeze(outliers_from_latent(cm_predict, eps=eps)) 
            n_outlier = len(outliers)
        except Exception as e:
            print(e)
            print("No outliers found")
            pass

        print(f'dimension = {model_dim}, eps = {eps}, number of outlier found: {n_outlier}')
        # get up to 1500 outliers 
        # if n_outlier > 1500:

        if n_outlier > outlier_max: 
            eps = eps + 0.09*random.random()
        elif n_outlier < outlier_min:
            eps = max(0.01, eps - 0.09*random.random())
        else: 
            eps_record[model_best] = eps 
            outlier_list.append(outliers) 
            break 
        outlier_count -= 1

    mytimer.mytime_label("outlier_eps", -1)
    return eps

def write_outliers(outlier_list, client):
    outlier_list_uni, outlier_count = np.unique(np.hstack(outlier_list), return_counts=True) 
    print(f"len(outlier_list_uni) = {len(outlier_list_uni)}, outlier_count = {outlier_count}")
    outliers_pdb_path = os.path.abspath('./outlier_pdbs') 
    make_dir_p(outliers_pdb_path) 
    print( 'Writing outliers in %s' % outliers_pdb_path  )

    new_outliers_list = [] 

    mytimer.mytime_label("outlier_pdb", 1)

    futures = []
    for outlier in outlier_list_uni:
        futures.append(client.submit(write_pdb, cvae_input[1][outlier], 
                                     cvae_input[2][outlier], 
                                     cvae_input[3][outlier], pdb_file, outliers_pdb_path))
    wait(futures)
    # garbage collect futures
    while(len(futures) > 0):
        del futures[0]

    mytimer.mytime_label("outlier_pdb", -1)

    for outlier in outlier_list_uni:
        myframe = cvae_input[1][outlier]
        myframe_v = cvae_input[3][outlier]
        hash = cvae_input[2][outlier]
        outlier_pdb_file = f'{outliers_pdb_path}/{hash}.pdb'
        outlier_v_file = f'{outliers_pdb_path}/{hash}.npy'
        new_outliers_list.append(outlier_pdb_file) 

    return new_outliers_list


def compute_rmsd(ref_pdb_file, restart_pdbs):
    mytimer.mytime_label("outlier_rmsd", 1)
    print("ref_pdf_file = ", ref_pdb_file)
    print("restart_pdbs[0] = ", restart_pdbs[0])
    print("len(restart_pdbs) = ", len(restart_pdbs))
    while(True):
        try:
            outlier_traj = mda.Universe(restart_pdbs[0], restart_pdbs) 
            break
        except Exception as e:
            print("Crashing while computing RMSD")
            print(e)
            time.sleep(3)
    ref_traj = mda.Universe(ref_pdb_file) 
    R = RMSD(outlier_traj, ref_traj, select='protein and name CA') 
    R.run()    
    restart_pdbs1 = [(rmsd, pdb) for rmsd, pdb in sorted(zip(R.rmsd[:,2], restart_pdbs))] 
    mytimer.mytime_label("outlier_rmsd", -1)
    return restart_pdbs1

def from_tmp():
    mytimer.mytime_label("outlier_lock", 1)

    dbfn1 = os.path.abspath('../') + "/OutlierDB.pickle"
    subprocess.getstatusoutput(f"touch {dbfn1}")

    mylock1 = LockFile(dbfn1)
    mylock1.acquire()

    print(subprocess.getstatusoutput("mv ../*.pickle ../pickle.tmp"))
    print(subprocess.getstatusoutput("mv ../outlier_pdbs ../tmp_outlier_pdbs"))
    print(subprocess.getstatusoutput("mv outlier_pdbs ../"))
    print(subprocess.getstatusoutput("mv *.pickle ../"))
    print(subprocess.getstatusoutput("mv *.json ../"))
    print(subprocess.getstatusoutput("rm -rf ../pickle.tmp ../tmp_outlier_pdbs"))

    mylock1.release()
    mytimer.mytime_label("outlier_lock", -1)
    return

def write_db(restart_pdb, restart_pdbs1):
    outlier_db_fn = os.path.abspath('OutlierDB.pickle')
    db = OutlierDB(os.path.abspath('../'), restart_pdbs1)
    with open(outlier_db_fn, 'wb') as f:
        pickle.dump(db, f)    
    restart_points = restart_pdbs
    restart_points_filepath = os.path.abspath('./restart_points.json') 
    restart_points = list(map(lambda x: x.replace("/tmp",""), restart_points))
    with open(restart_points_filepath, 'w') as restart_file: 
        json.dump(restart_points, restart_file) 
    return db

if __name__ == '__main__':
    pprint(dict(os.environ))
    DEBUG = 1 
    args = parse_arguments()

    pdb_file = os.path.abspath(args.pdb) 
    ref_pdb_file = os.path.abspath(args.ref) 

    waiting_time1 = 60
    min_step_increment = 100
    waiting_time2 = 60
    current_step = 0
    waiting_time3 = 60

    outlier_max = 4000
    outlier_min = 3000
    outlier_count = 120

    lastN = 8_000

    bp = os.path.join(args.md, '*.bp')
    print("bp=",bp)

    print("hostname = ", subprocess.getstatusoutput("hostname")[1])
    adios_files_list = wait_for_bp(args.md, args.n, waiting_time1)
    steps = wait_for_bp_nonempty(adios_files_list, waiting_time1)
    steps = wait_for_enough_data(adios_files_list, steps, waiting_time2, current_step, min_step_increment)

    client = Client(processes=True, n_workers=39, local_directory='/tmp')
    print("Client = ", client)

    mystreams = STREAMS(adios_files_list, lastN=lastN)
    j = 0

    while(True):
        print(f"outlier iteration {j}")
        mytimer.mytime_label("outlier_search", 1)

        top_dir, tmp_dir = run_dir()
        model_best = wait_for_model(args.cvae, waiting_time3)
        cvae_input = outlier_read(mystreams)
        eps, eps_record, model_dim, eps_record_filepath = eps_init(model_best)
        cm_predict = predict(model_best, cvae_input, model_dim)
        outlier_list = []
        cluster(cm_predict, outlier_list, eps_record, eps, model_best,
                outlier_count, outlier_max, outlier_min)

        with open(eps_record_filepath, 'w') as eps_file: 
            json.dump(eps_record, eps_file)
 
        restart_pdbs = write_outliers(outlier_list, client)

        if(len(restart_pdbs) == 0):
            print("No outliers found")
            os.chdir(top_dir)
            mytimer.mytime_label("outlier_search", -1)
            j += 1
            continue

        restart_pdbs1 = compute_rmsd(ref_pdb_file, restart_pdbs)
        db = write_db(restart_pdbs, restart_pdbs1)
        from_tmp()

        os.chdir(top_dir)
        mytimer.mytime_label("outlier_search", -1)
        j += 1


