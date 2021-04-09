import simtk.openmm.app as app
import simtk.openmm as omm
import simtk.unit as u 
import os
import mytimer
import numpy as np 
import h5py 
import sys

import adios2
import hashlib

from MDAnalysis.analysis import distances
from myconvert import *

step = 0

class ContactMapReporter(object):
    def __init__(self, file, reportInterval):
        #self._file = h5py.File(file, 'w', libver='latest')
        #self._file.swmr_mode = True
        #self._out = self._file.create_dataset('contact_maps', shape=(2,0), maxshape=(None, None))
        self._reportInterval = reportInterval
        self._adios_stream = adios2.open(name="cms_positions.bp", mode="w", config_file="adios.xml", io_in_config_file=os.path.basename(os.getcwd()))
        print("ContactMapReporter constructor"); sys.stdout.flush()
    def __del__(self):
        #self._file.close()
        self._adios_stream.close()
        print("ContactMapReporter destructor"); sys.stdout.flush()
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, None)

    def report(self, simulation, state):
        mytimer.mytime_label("simulation_report",start=1)
        global step
        #print(f"type(simulation) = {type(simulation)}, dir(simulation) = {dir(simulation)}")
        stateA = simulation.context.getState(getPositions=True, getVelocities=True)
        ca_indices = []
        for atom in simulation.topology.atoms():
            if atom.name == 'CA':
                ca_indices.append(atom.index)
        positions = np.array(state.getPositions().value_in_unit(u.angstrom))
        velocities = stateA.getVelocities(asNumpy=True)
        #print(f"type(velocities) = {type(velocities)}, velocities.shape = {velocities.shape}, dir(velocities) = {dir(velocities)}")
        #print(f"type(velocities._value) = {type(velocities._value)}")
        #print(f"type(velocities.value_in_unit) = {type(velocities.value_in_unit)}")

        sys.stdout.flush()

        m = hashlib.md5()
        m.update(positions.tostring())
        md5 = m.hexdigest()

        '''
        print(f'type(md5) = {type(md5)}')
        print(f'md5 = {md5}')
        '''

        md5 = hash2intarray(md5)


        time = int(np.round(state.getTime().value_in_unit(u.picosecond)))
        positions_ca = positions[ca_indices].astype(np.float32)
        distance_matrix = distances.self_distance_array(positions_ca)
        contact_map = np.asarray((distance_matrix < 8.0), dtype=np.int32)
        stepA = np.array([step], dtype=np.int32)

        # print(f"reporting: time = {time}, step = {step}")


        '''
        print(f"md5.dtype = {md5.dtype}")
        print(f"stepA.dtype = {stepA.dtype}")
        print(f"positions.dtype = {positions.dtype}")
        print(f"contact_map.dtype = {contact_map.dtype}")
        print(f"positions_ca.shape = {positions_ca.shape}")
        print(f"distance_matrix.shape = {distance_matrix.shape}")
        print(f"contact_map.shape = {contact_map.shape}")
        print(f"positions.shape = {positions.shape}")
        sys.stdout.flush()
        '''

        self._adios_stream.write("md5", md5, list(md5.shape), [0]*len(md5.shape), list(md5.shape))
        self._adios_stream.write("step", stepA, list(stepA.shape), [0]*len(stepA.shape), list(stepA.shape))
        self._adios_stream.write("positions", positions, list(positions.shape), [0]*len(positions.shape), list(positions.shape))
        self._adios_stream.write("velocities", velocities._value, list(velocities.shape), [0]*len(velocities.shape), list(velocities.shape))
        self._adios_stream.write("contact_map", contact_map, list(contact_map.shape), [0]*len(contact_map.shape), list(contact_map.shape), end_step=True)
        '''
        new_shape = (len(contact_map), self._out.shape[1] + 1) 
        print("ContactMapReporter: In report, time = ", time); sys.stdout.flush()
        self._out.resize(new_shape)
        self._out[:, new_shape[1]-1] =contact_map
        self._file.flush()
        '''
        step += 1
        mytimer.mytime_label("simulation_report",start=-1)
