import adios2
import time
import random
import numpy as np
from myutils import *
import os.path
import sys
import subprocess

class Simulation:
    def __init__(self, dir, adios_xml, aggregator_dir):
        self.dir = dir
        self.adios_xml = adios_xml
        self.aggregator_dir = aggregator_dir
        subprocess.getstatusoutput(f'mkdir -p {dir}')
        subprocess.getstatusoutput(f'cp {adios_xml} {dir}/adios.xml')
        new_dir = dir.replace("all","new")
        print(f'ln -s {dir} {new_dir}')
        subprocess.getstatusoutput(f'ln -s {dir} {new_dir}')
        self.log = open(f'{dir}/simulation.log', 'w')
        self.log.write("Start\n")
        self.log.write(f'{get_now()}\n')
        self.step = 0
        self._adios_stream = adios2.open(name=f"{dir}/SimulationOutput.bp", mode="w", config_file=f"{dir}/adios.xml", io_in_config_file="SimulationOutput") 
        self.stop = False
    def __del__(self):
        self._adios_stream.close()
        self.log.close()
    def produce(self):
        self.data = np.random.rand(3,2)
    def iterate(self):
        self.log.write(f'In iterate: step = {self.step}\n')
        self.produce()
        self.log.write(f'data = {self.data}\n')
        self.step += 1
        self._adios_stream.write("MyData", self.data, list(self.data.shape), [0,0], list(self.data.shape), end_step=True)
        self.qstop()
    def qstop(self):
        if(os.path.exists(f"{self.dir}/stop.simulation") or os.path.exists(f"{self.aggregator_dir}/stop.aggregator")):
           self.log.write("Received kill signal\n")
           self.stop = True
    def set_stop(self):
        self.stop = True
    def run(self):
        while(not self.stop):
            self.log.write("="*30 + "\n")
            delay = random.uniform(0,3)
            self.log.write(f'Sleeping for {delay:.1f}\n')
            time.sleep(delay)
            self.iterate()

if(__name__ == '__main__'):
    try:
        sim_dir = sys.argv[1]        
        adios_xml = sys.argv[2]
        task = sys.argv[3]
        aggregator_dir = sys.argv[4]
    except:
        print("Usage: python simulation.py dir adios_xml task")
        sys.exit(1)
    counter = 0
    while(not os.path.exists(f"{aggregator_dir}/stop.aggregator")):
        dir = f'{sim_dir}/{task}_{counter}'
        s = Simulation(dir, adios_xml, aggregator_dir)
        s.run()
        del s
        counter += 1



    
