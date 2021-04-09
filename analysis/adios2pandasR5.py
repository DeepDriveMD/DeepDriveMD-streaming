#!/usr/bin/env python

from radical.entk import Pipeline, Stage, Task, AppManager
import os

# ------------------------------------------------------------------------------
# Set default verbosity

if os.environ.get('RADICAL_ENTK_VERBOSE') is None:
    os.environ['RADICAL_ENTK_REPORT'] = 'True'

# Description of how the RabbitMQ process is accessible
# No need to change/set any variables if you installed RabbitMQ has a system
# process. If you are running RabbitMQ under a docker container or another
# VM, set "RMQ_HOSTNAME" and "RMQ_PORT" in the session where you are running
# this script.
hostname = os.environ.get('RMQ_HOSTNAME', 'localhost')
port = int(os.environ.get('RMQ_PORT', 5672))
username = os.environ.get('RMQ_USERNAME')
password = os.environ.get('RMQ_PASSWORD')


queue = 'batch'
proj_id = 'csc299'


if __name__ == '__main__':

    # Create a Pipeline object
    p = Pipeline()

    # Create a Stage object
    s = Stage()

    # Create a Task object
    t = Task()
    t.pre_exec = ['. /ccs/home/iyakushin/etc/openmm5.sh']
    t.pre_exec += ['cd /gpfs/alpine/scratch/iyakushin/csc299/Test3/entk_cvae_md_devel3/entk_cvae_md/analysis']
    t.name = 'test'
    t.executable = ['/ccs/home/iyakushin/.conda/envs/OpenMM5/bin/python']
    t.arguments = ['/gpfs/alpine/scratch/iyakushin/csc299/Test3/entk_cvae_md_devel3/entk_cvae_md/analysis/adios2pandasR4.py', 
                   '/gpfs/alpine/scratch/iyakushin/csc299/Test3/entk_cvae_md_devel3/entk_cvae_md/aggregate/aggregator0.bp', 42]


    t.cpu_reqs = {'processes': 42,
                  'process_type': None,
                  'threads_per_process': 4,
                  'thread_type': 'OpenMP'
            }
    t.gpu_reqs = {'processes': 1,
                  'process_type': None,
                  'threads_per_process': 1,
                  'thread_type': 'CUDA'
              }

    # Add Task to the Stage
    s.add_tasks(t)

    # Add Stage to the Pipeline
    p.add_stages(s)

    res_dict = {
            'resource': 'ornl.summit',
            'queue'   : queue,
            'schema'  : 'local',
            'walltime': 60,
            'cpus'    : 42,
            'gpus'    : 1,
            'project' : proj_id
    }


    # Create Application Manager
    # appman = AppManager()
    appman = AppManager(hostname=os.environ.get('RMQ_HOSTNAME'), port=int(os.environ.get('RMQ_PORT')), 
                        username=os.environ.get('RMQ_USERNAME'), password=os.environ.get('RMQ_PASSWORD'))

    appman.resource_desc = res_dict

    appman.workflow = set([p])

    # Run the Application Manager
    appman.run()



"""

import numpy as np
import adios2
import MDAnalysis as mda
from  MDAnalysis.analysis.rms import RMSD
import sys
import pandas as pd
from multiprocessing import Pool
import glob
import os

fn = sys.argv[1]
ps = int(sys.argv[2])

def f(position):
    outlier_traj = mda.Universe(init_pdb, position)
    ref_traj = mda.Universe(ref_pdb_file)
    R = RMSD(outlier_traj, ref_traj, select = 'protein and name CA')
    R.run()
    return R.rmsd[:,2][0]

ref_pdb_file ='../MD_exps/fs-pep/pdb/1FME-0.pdb'

init_pdb = '../MD_exps/fs-pep/pdb/1FME.pdb'

pf = pd.DataFrame(columns=["fstep", "step", "R"])

with adios2.open(fn, "r") as fr:
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
    
    name = 'step'
    steps = fr.read(name, [],[], 0, n)
    print(type(steps))
    print(steps.shape)
    print(steps)
    sys.stdout.flush()


with Pool(processes=ps) as pool:
    Rs = pool.map(f, positions)

pf['fstep'] = list(np.arange(len(steps)))
pf['step'] = steps
pf['R'] = Rs

fn1 = os.path.basename(fn).replace(".bp",".csv")

pf.to_csv(fn1)



            t1.pre_exec = ['. /ccs/home/iyakushin/etc/openmm5.sh']
            t1.pre_exec += ['export PYTHONPATH=%s/MD_exps:$PYTHONPATH' % base_path] 
            t1.pre_exec += ['cd %s' % md_path]
            t1.executable = ['%s/bin/python' % conda_path]  # run_openmm.py
            # t1.arguments = ['-m', 'cProfile', '-o', 'run_openmm.cprofile', '%s/run_openmm.py' % md_path] 
            t1.arguments = ['%s/run_openmm.py' % md_path] 
            if top_file: 
                t1.arguments += ['--topol', top_file]

            t1.arguments += ['-d', base_path]
            t1.arguments += ['-t', i]

            # pick initial point of simulation 
            if initial_MD or i >= len(outlier_list): 
                t1.arguments += ['--pdb_file', pdb_file] 
                t1.arguments += ['-i']
            elif outlier_list[i].endswith('pdb'): 
                t1.arguments += ['--pdb_file', outlier_list[i]] 
            elif outlier_list[i].endswith('chk'): 
                t1.arguments += ['--pdb_file', pdb_file, 
                        '-c', outlier_list[i]] 
            if initial_MD: 
                t1.arguments += ['--length', LEN_initial] 
            else: 
                t1.arguments += ['--length', LEN_iter]

            # assign hardware the task 
            t1.cpu_reqs = {'processes': 1,
                           'process_type': None,
                              'threads_per_process': 4,
                              'thread_type': 'OpenMP'
                              }



"""
