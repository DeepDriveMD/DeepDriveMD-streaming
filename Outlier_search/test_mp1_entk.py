from radical.entk import Pipeline, Stage, Task, AppManager
import os


envfn = '/g/g15/yakushin/etc/powerai.sh'
PYTHON = '/g/g15/yakushin/.conda/envs/powerai/bin/python'
dir = '/p/gpfs1/yakushin/DDMD/22/entk_cvae_md/Outlier_search'

def generate_pipeline():
    p = Pipeline()
    p.name = 'MP_pipeline'

    s = Stage()
    s.name = 'MP_stage'

    p.add_stages(s)

    t = Task()
    t.name = 'MP_task'
    t.pre_exec = ['source ' + envfn]
    t.pre_exec = ['cd ' + dir]
    t.executable = [PYTHON]
    t.arguments = [dir + '/test_mp1.py'] 
    
    t.cpu_reqs = {'processes': 1,
                   'process_type': None,
                   'threads_per_process': 156,
                   'thread_type': 'OpenMP'
    }
    t.gpu_reqs = {'processes': 1,
                   'process_type': None,
                   'threads_per_process': 1,
                   'thread_type': 'CUDA'
    }
    
    # Add the learn task to the learning stage
    s.add_tasks(t)

    return [p]

queue = 'pbatch'
proj_id = 'cv19-a01'
hrs_wt = 1


res_dict = {
    'resource': 'llnl.lassen',
    'queue'   : queue,
    'schema'  : 'local',
    'walltime': 60 * hrs_wt,
    'cpus'    : 159,
    'gpus'    : 1,
    'project' : proj_id
}
    
    
appman = AppManager(hostname=os.environ.get('RMQ_HOSTNAME'), port=int(os.environ.get('RMQ_PORT')), 
                    username=os.environ.get('RMQ_USERNAME'), password=os.environ.get('RMQ_PASSWORD'))
    
appman.resource_desc = res_dict    
pipelines = generate_pipeline()
appman.workflow = pipelines
appman.run()
