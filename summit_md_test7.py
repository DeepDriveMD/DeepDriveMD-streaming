import os, json, time 
from radical.entk import Pipeline, Stage, Task, AppManager
import subprocess
import math

# ------------------------------------------------------------------------------
# Set default verbosity

if os.environ.get('RADICAL_ENTK_VERBOSE') is None:
    os.environ['RADICAL_ENTK_REPORT'] = 'True'

# Assumptions:
# - # of MD steps: 2
# - Each MD step runtime: 15 minutes
# - Summit's scheduling policy [1]
#
# Resource rquest:
# - 4 <= nodes with 2h walltime.
#
# Workflow [2]
#
# [1] https://www.olcf.ornl.gov/for-users/system-user-guides/summit/summit-user-guide/scheduling-policy
# [2] https://docs.google.com/document/d/1XFgg4rlh7Y2nckH0fkiZTxfauadZn_zSn3sh51kNyKE/
#
'''
export RMQ_HOSTNAME=two.radical-project.org 
export RMQ_PORT=33235 
export RADICAL_PILOT_DBURL=mongodb://hyperrct:h1p3rrc7@two.radical-project.org:27017/hyperrct 
export RADICAL_PROFILE=True
export RADICAL_PILOT_PROFILE=True
export RADICAL_ENTK_PROFILE=True
'''
#

home_path = os.getenv('HOME')
#env_fn = home_path + '/etc/openmm_14.sh
env_fn = "/usr/workspace/cv_ddmd/software1/etc/powerai.sh"
print(env_fn)
# env_fn1 = home_path + '/etc/powerai.sh'
env_fn1 = "/usr/workspace/cv_ddmd/software1/etc/powerai.sh"
print(env_fn1)

base_path = os.path.abspath('.') # '/gpfs/alpine/proj-shared/bip179/entk/hyperspace/microscope/experiments/'
misc_path = f"{base_path}/misc"
#conda_path = home_path + "/.conda/envs/DDMD5/"
#conda_path1 = home_path + "/.conda/envs/powerai/"

conda_path = "/usr/workspace/cv_ddmd/conda1/powerai/"
conda_path1 = "/usr/workspace/cv_ddmd/conda1/powerai/"


print(f"base_path = {base_path}")
print(f"conda_path = {conda_path}")

md_path = os.path.join(base_path, 'MD_exps/fs-pep') 
#agg_path = os.path.join(base_path, 'MD_to_CVAE') 
agg_path = os.path.join(base_path, 'aggregate') 
cvae_path = os.path.join(base_path, 'CVAE_exps') 
outlier_path = os.path.join(base_path, 'Outlier_search') 

agg_md_path = os.path.join(base_path, 'aggregate') 

barrier_fn = f'{base_path}/b1.barrier'


# pdb_file = os.path.join(md_path, 'pdb/100-fs-peptide-400K.pdb') 
pdb_file = os.path.join(md_path, 'pdb/1FME.pdb') 
top_file = None 
# ref_pdb_file = os.path.join(md_path, 'pdb/fs-peptide.pdb')
ref_pdb_file = os.path.join(md_path, 'pdb/1FME-0.pdb')

'''
N_jobs_MD = 12
'''

N_jobs_MD = 120
N_jobs_ML = 1
dim_offset_ML = 10
N_aggregators = 10
hrs_wt = 12
LEN_initial = 10 # 100
#LEN_initial = 10

'''
### debug start
N_jobs_MD=2
N_jobs_ML=1
dim_offset_ML = 10
N_aggregators=1
hrs_wt = 1
LEN_initial = 10
### debug end
'''

queue = 'pbatch'

proj_id = 'cv19-a01'

CUR_STAGE=0
MAX_STAGE= 2
RETRAIN_FREQ = 5



LEN_iter = 10 

'''
LEN_initial = 10
LEN_iter = 10
'''


def generate_training_pipeline():
    """
    Function to generate the CVAE_MD pipeline
    """

    def generate_clean_stage():
        s = Stage()
        s.name='clean'
        t = Task()
        t.pre_exec = [f'. {env_fn}']
        t.pre_exec += ['cd %s' % base_path]
        t.executable = ['make']
        t.arguments = ['clean']
        t.post_exec = ['touch clean.done']
        t.cpu_reqs = {'processes': 1,
                      'process_type': None,
                      'threads_per_process': 4,
                      'thread_type': 'OpenMP'
        }

        s.add_tasks(t)
        return s

    def generate_wait4clean_stage():
        s = Stage()
        s.name='wait4clean'
        t = Task()
        t.pre_exec = [f'. {env_fn}']
        t.pre_exec += ['cd %s' % base_path]
        t.executable = ['%s/bin/python' % conda_path]
        t.arguments = ['wait4clean.py']
        t.arguments += ['clean.done']
        t.cpu_reqs = {'processes': 1,
                      'process_type': None,
                      'threads_per_process': 4,
                      'thread_type': 'OpenMP'
        }
        s.add_tasks(t)
        return s
        

    def generate_MD_stage(num_MD=1): 
        """
        Function to generate MD stage. 
        """
        s1 = Stage()
        s1.name = 'MD'
        initial_MD = True 
        outlier_filepath = '%s/restart_points.json' % outlier_path

        if os.path.exists(outlier_filepath): 
            initial_MD = False 
            outlier_file = open(outlier_filepath, 'r') 
            outlier_list = json.load(outlier_file) 
            print(f"Exists {outlier_filepath}")
            outlier_file.close() 

        # MD tasks
        time_stamp = int(time.time())
        for i in range(num_MD):
            t1 = Task()
            t1.pre_exec = [f'. {env_fn}']
            t1.pre_exec += [f'export PYTHONPATH={base_path}/MD_exps:{misc_path}:$PYTHONPATH'] 
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
            t1.gpu_reqs = {'processes': 1,
                           'process_type': None,
                              'threads_per_process': 1,
                              'thread_type': 'CUDA'
                             }
                              
            # Add the MD task to the simulating stage
            s1.add_tasks(t1)
            # s1.post_exec = post_MD
        return s1 

    def post_MD():
        subprocess.getstatusoutput(f'touch {barrier_fn}')


    def generate_aggregate_md_stage(start_unit, end_unit, output_postfix):
        s2 = Stage()
        s2.name = 'aggregating_md'

        # Aggregation task
        t2 = Task()
        t2.pre_exec = [] 
        t2.pre_exec += [f'. {env_fn}']
        t2.pre_exec += [f'export PYTHONPATH={misc_path}:$PYTHONPATH']
        t2.pre_exec += ['cd %s' % agg_md_path]
        t2.executable = ['%s/bin/python' % conda_path]  # MD_to_CVAE.py
        # t2.arguments = ['-m', 'cProfile', '-o', 'aggregate.cprofile', '%s/aggregator.py' % agg_md_path, base_path, base_path, 10000]
        t2.arguments = ['%s/aggregator.py' % agg_md_path, '-c', base_path, '-r', base_path, '-m', 100000, '-s', start_unit, '-e', end_unit, '-o', output_postfix]
        # t2.post_exec = [f'rm -f {base_path}/clean.done']

        t2.cpu_reqs = {'processes': 1,
                       'process_type': None,
                       'threads_per_process': 4,
                       'thread_type': 'OpenMP'
        }
        s2.add_tasks(t2)
        return s2         


    def generate_aggregating_stage(): 
        """ 
        Function to concatenate the MD trajectory (h5 contact map) 
        """ 
        s2 = Stage()
        s2.name = 'aggregating'

        # Aggregation task
        t2 = Task()
        # https://github.com/radical-collaboration/hyperspace/blob/MD/microscope/experiments/MD_to_CVAE/MD_to_CVAE.py
        t2.pre_exec = [] 
        t2.pre_exec += [f'. {env_fn}']
        t2.pre_exec += ['cd %s' % agg_path]
        t2.executable = ['%s/bin/python' % conda_path]  # MD_to_CVAE.py
        t2.arguments = ['%s/MD_to_CVAE.py' % agg_path, 
                '--sim_path', md_path]

        # Add the aggregation task to the aggreagating stage
        s2.add_tasks(t2)
        return s2 


    def generate_ML_stage(num_ML=1): 
        """
        Function to generate the learning stage
        """
        s3 = Stage()
        s3.name = 'learning'

        # learn task
        time_stamp = int(time.time())
        for i in range(num_ML): 
            t3 = Task()
            # https://github.com/radical-collaboration/hyperspace/blob/MD/microscope/experiments/CVAE_exps/train_cvae.py
            t3.pre_exec = []
            t3.pre_exec += [f'. {env_fn1}']
            t3.pre_exec += ['export PYTHONPATH=%s/CVAE_exps:%s/misc:$PYTHONPATH' % (base_path, base_path)]
            t3.pre_exec += ['cd %s' % cvae_path]
            dim = i + dim_offset_ML 
            cvae_dir = 'cvae_runs_%.2d_%d' % (dim, time_stamp+i) 
            t3.pre_exec += ['mkdir -p {0} && cd {0}'.format(cvae_dir)]
            t3.executable = ['%s/bin/python' % conda_path1]  # train_cvae.py
            t3.arguments = ['%s/train_cvae.py' % cvae_path, 
                    '--bp_dir', agg_path, '--nbp_files', N_aggregators,  
                    '--dim', dim] 
            
            t3.cpu_reqs = {'processes': 1,
                           'process_type': None,
                    'threads_per_process': 4,
                    'thread_type': 'OpenMP'
                    }
            t3.gpu_reqs = {'processes': 1,
                           'process_type': None,
                    'threads_per_process': 1,
                    'thread_type': 'CUDA'
                    }
        
            # Add the learn task to the learning stage
            s3.add_tasks(t3)

        return s3 


    def wait_4_b1():
        s4w = Stage()
        s4w.name = 'wait for b1'

        t4w = Task()
        t4w.pre_exec = [] 
        t4w.pre_exec += [f'. {env_fn}']
        t4w.pre_exec += ['export PYTHONPATH=%s/CVAE_exps:$PYTHONPATH' % base_path] 
        t4w.executable = ['%s/bin/python' % conda_path] 
        t4w.arguments = ['%s/wait_4_b1.py' % base_path, barrier_fn ]

        t4w.cpu_reqs = {'processes': 1,
                        'process_type': None,
                        'threads_per_process': 4,
                        'thread_type': 'OpenMP'
                }
        t4w.gpu_reqs = {'processes': 0,
                           'process_type': None,
                'threads_per_process': 0,
                'thread_type': 'CUDA'
                }
        s4w.add_tasks(t4w) 
        return s4w

    def generate_interfacing_stage(simplify=False): 
        s4 = Stage()
        s4.name = 'scanning'

        # Scaning for outliers and prepare the next stage of MDs 
        t4 = Task() 
        t4.pre_exec = [] 
        t4.pre_exec += [f'. {env_fn1}']
        t4.pre_exec += ['export PYTHONPATH=%s/CVAE_exps:%s/misc:$PYTHONPATH' % (base_path, base_path)] 
        t4.pre_exec += ['cd %s/Outlier_search' % base_path]
        # t4.pre_exec += ['export OMP_NUM_THREADS=32']
        # t4.pre_exec += ['export POWERAI_SAVE_OMP_NUM_THREADS=32']
        t4.pre_exec += ['export OMP_PROC_BIND=false']
        t4.executable = ['%s/bin/python' % conda_path1] 
        t4.arguments = ['outlier_locator.py', 
                        '--md',  agg_path,
                        '--nbp_files', N_aggregators,
                        '--cvae', cvae_path, 
                        '--pdb', pdb_file, 
                        '--ref', ref_pdb_file]

        t4.cpu_reqs = {'processes': 1,
                       'process_type': None,
                       'threads_per_process': 39,
                       'thread_type': 'OpenMP'
                   }
        t4.gpu_reqs = {'processes': 1,
                       'process_type': None,
                       'threads_per_process': 1,
                       'thread_type': 'CUDA'
                   }
        s4.add_tasks(t4) 
        #if(not simplify):
        #    s4.post_exec = func_condition 
        
        return s4


    def func_condition(): 
        global CUR_STAGE, MAX_STAGE 
        if CUR_STAGE < MAX_STAGE: 
            func_on_true()
        else:
            func_on_false()

    def func_on_true(): 
        global CUR_STAGE, MAX_STAGE
        print('finishing stage %d of %d' % (CUR_STAGE, MAX_STAGE)) 
        
        # --------------------------
        # MD stage
        s1 = generate_MD_stage(num_MD=N_jobs_MD)
        # Add simulating stage to the training pipeline
        p.add_stages(s1)

        if CUR_STAGE % RETRAIN_FREQ == 0: 
            # --------------------------
            # Aggregate stage
            s2 = generate_aggregating_stage() 
            # Add the aggregating stage to the training pipeline
            p.add_stages(s2)

            # --------------------------
            # Learning stage
            s3 = generate_ML_stage(num_ML=N_jobs_ML) 
            # Add the learning stage to the pipeline
            p.add_stages(s3)

        # --------------------------
        # Outlier identification stage
        s4 = generate_interfacing_stage() 
        p.add_stages(s4) 
        
        CUR_STAGE += 1

    def func_on_false(): 
        print('Done') 


    def test4(): 
        s4 = generate_interfacing_stage(simplify=True)
        p.add_stages(s4) 


    def test5():
        s1 = generate_MD_stage(num_MD=N_jobs_MD)
        p_md.add_stages(s1)

        '''
        s2 = generate_aggregating_stage()
        p_train.add_stages(s2)
        s3 = generate_ML_stage(num_ML=N_jobs_ML)
        p_train.add_stages(s3)

        '''
        s4w = wait_4_b1()
        s4 = generate_interfacing_stage(True)
        p_outliers.add_stages(s4w)
        p_outliers.add_stages(s4)

    def test6():
        s1 = generate_MD_stage(num_MD=N_jobs_MD)
        p_md.add_stages(s1)

        s2 = generate_aggregate_md_stage()
        p_outliers.add_stages(s2)


    def test7():
        s1 = generate_MD_stage(num_MD=N_jobs_MD)
        p_md.add_stages(s1)

        s2 = generate_aggregate_md_stage()
        p_outliers.add_stages(s2)
        s3 = generate_ML_stage(num_ML=N_jobs_ML)
        p_outliers.add_stages(s3)

    def test8():
        s1 = generate_MD_stage(num_MD=N_jobs_MD)
        p_md.add_stages(s1)

        s2 = generate_aggregate_md_stage()
        p_outliers.add_stages(s2)
        s3 = generate_ML_stage(num_ML=N_jobs_ML)
        p_outliers.add_stages(s3)
        s4 = generate_interfacing_stage(simplify=True)
        p_outliers.add_stages(s4)


    def test9():
        #p_md.add_stages(generate_clean_stage())
        


        s1 = generate_MD_stage(num_MD=N_jobs_MD)
        p_md.add_stages(s1)

        step = math.ceil(N_jobs_MD/N_aggregators)

        for k in range(N_aggregators):
            s2 = generate_aggregate_md_stage(k*step, min((k+1)*step, N_jobs_MD) - 1, k)
            #p_aggregate[k].add_stages(generate_wait4clean_stage())
            p_aggregate[k].add_stages(s2)

        #p_ml.add_stages(generate_wait4clean_stage())
        s3 = generate_ML_stage(num_ML=N_jobs_ML)
        p_ml.add_stages(s3)

        #p_outliers.add_stages(generate_wait4clean_stage())
        s4 = generate_interfacing_stage(simplify=True)
        p_outliers.add_stages(s4)

    global CUR_STAGE

    pipelines = []

    p_md = Pipeline()
    p_md.name = 'MD'

    pipelines.append(p_md)

    p_aggregate = []
    for k in range(N_aggregators):
        p_aggregate.append(Pipeline())
        p_aggregate[k].name = 'aggregate_%d' % k
        pipelines.append(p_aggregate[k])

    p_ml = Pipeline()
    p_ml.name = 'ML'
    pipelines.append(p_ml)

    p_outliers = Pipeline()
    p_outliers.name = 'outliers'
    pipelines.append(p_outliers)

    test9()

    #    func_on_true()
    
    return pipelines


if __name__ == '__main__':

    # Create a dictionary to describe four mandatory keys:
    # resource, walltime, cores and project
    # resource is 'local.localhost' to execute locally
    '''
    res_dict = {
            'resource': 'ornl.summit',
            'queue'   : queue,
            'schema'  : 'local',
            'walltime': 60 * hrs_wt,
            'cpus'    : N_jobs_MD * 7 + 3,
            'gpus'    : N_jobs_MD + 1,#6*2 ,
            'project' : proj_id,
      }
    '''

    res_dict = {
        'resource': 'llnl.lassen',
        'queue'   : queue,
        'schema'  : 'local',
        'walltime': 60 * hrs_wt,
        'cpus'    : 4*(N_jobs_MD + N_jobs_ML + N_aggregators + 39) + 160*4,
        'gpus'    : N_jobs_MD + N_jobs_ML  + 4,
        'project' : proj_id
    }
    
    
    # Create Application Manager
    # appman = AppManager()
    appman = AppManager(hostname=os.environ.get('RMQ_HOSTNAME'), port=int(os.environ.get('RMQ_PORT')), 
                        username=os.environ.get('RMQ_USERNAME'), password=os.environ.get('RMQ_PASSWORD'))
    
    appman.resource_desc = res_dict
    
    pipelines = generate_training_pipeline()
    
    appman.workflow = pipelines
    
    # Run the Application Manager
    appman.run()
