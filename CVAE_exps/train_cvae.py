import os, sys, errno
import argparse 
from cvae.CVAE import run_cvae  
import numpy as np 
import time
import subprocess
import mytimer
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--bp_dir", dest="f", default='', help="Aggregate directory with bp files")
parser.add_argument("-n", "--nbp_files", dest="n", default=1, type=int, help="Number of bp files")
# parser.add_argument("-o", help="output: cvae weight file. (Keras cannot load model directly, will check again...)")
parser.add_argument("-d", "--dim", default=3, help="Number of dimensions in latent space")
parser.add_argument("-gpu", default=0, help="gpu_id")

args = parser.parse_args()

#cvae_input = glob.glob(args.f + '*.bp')

hyper_dim = int(args.dim) 
gpu_id = args.gpu

#h5_input = args.h

print("hostname = ", subprocess.getstatusoutput("hostname")[1])

min_step_increment = 5000
max_steps = 400000

#min_step_increment = 50

current_step = 0

i = 0

loss = 10000

while(True):
    print(f"training iteration {i}")

    while(True):
        cvae_input = mytimer.get_bp_files(args.f, args.n)
        if(cvae_input == None):
            print(f"Waiting for {cvae_input} to be produced")
            sys.stdout.flush()
            time.sleep(60)
        else:
            break

    print(f"cvae_input={cvae_input}")

    while(True):
        try:
            steps = mytimer.bp_steps(cvae_input)
            print(f"steps = {steps}")
            sys.stdout.flush()
            break
        except Exception as e:
            print(e)
            print(f"Waiting for {cvae_input} to become non-empty")
            sys.stdout.flush()
        time.sleep(60)


    while(steps < min_step_increment):
        print(f"Waiting for at least {min_step_increment}, currently steps = {steps}")
        sys.stdout.flush()
        time.sleep(60)
        steps = mytimer.bp_steps(cvae_input)


    print(f'steps={steps}')
    sys.stdout.flush()

    '''
    while(steps - current_step < min_step_increment):
        print(f"Waiting for enough time steps to accumulate since the last time: {steps} - {current_step} < {min_step_increment}")
        sys.stdout.flush()
        time.sleep(60)
        steps = int(subprocess.getstatusoutput(f"bpls {cvae_input}")[1].split('\n')[0].split('*')[0].split(" ")[-1])
    '''

    current_step = steps

    top_dir = os.getcwd()

    print(f"top_dir={top_dir}")

    tmp_dir = 'tmp'
    if(not os.path.exists(tmp_dir)):
        os.mkdir(tmp_dir)
    os.chdir(tmp_dir)

    mytimer.mytime_label("cvae_iteration",1)

    # mytimer.mytop()

    print(f"loss = {loss}")
    if(loss > 1000):
        '''
        try:
            cvae
        except:
            cvae = run_cvae(gpu_id, cvae_input, hyper_dim=hyper_dim, nsamples=max_steps, epochs=1)
        cvae.load('../../cvae_init/best.h5')
        '''
        cvae = run_cvae(gpu_id, cvae_input, hyper_dim=hyper_dim, nsamples=max_steps, epochs=20)        
    else:
        cvae = run_cvae(gpu_id, cvae_input, hyper_dim=hyper_dim, nsamples=max_steps, init_cvae = cvae, epochs=20)

    mytimer.mytime_label("cvae_iteration",-1)

    print(f"currend_dir = {os.getcwd()}"); 
    print(subprocess.getstatusoutput("ls -l"))
    sys.stdout.flush()

    if(not os.path.exists('best.h5')):
        os.chdir(top_dir)
        cvae = None
        i += 1
        continue


    cvae.load('best.h5')


    model_weight = r'cvae_weight.h5' 
    model_file = r'cvae_model.h5'
    model_json = r'cvae_model.json'
    loss_file = r'loss.npy' 

    cvae.model.save_weights(model_weight)
    cvae.save(model_file)
    f = open(model_json, "w")
    f.write(cvae.model.to_json())
    f.close()
    np.save(loss_file, cvae.history.val_losses) 

    loss = cvae.history.val_losses[-1]

    # del cvae

    #subprocess.getstatusoutput(f"mkdir -p ../archive/{i}")
    #subprocess.getstatusoutput(f"cp -a * ../archive/{i}")

    subprocess.getstatusoutput('mv * ../')
    os.chdir(top_dir)
    i += 1
