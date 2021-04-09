import simtk.openmm.app as app
import simtk.openmm as omm
import simtk.unit as u
import os
import subprocess
import mytimer
import sys
import pickle

from OutlierDB import *
from lockfile import LockFile

import parmed as pmd
import random
from MD_utils.openmm_reporter import ContactMapReporter


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
        # pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/100-fs-peptide-400K.pdb'
        pdb_file = f'{top_dir}/MD_exps/fs-pep/pdb/1FME.pdb'
    return pdb_file

def mytop():
    print("="*10 + " top " + "="*10)
    print(subprocess.getstatusoutput("top -U iyakushin -b -n 1")[1])
    print("="*25)

def openmm_simulate_charmm_nvt(top_file, pdb_file, check_point=None, GPU_index=0,  
        output_traj="output.dcd", output_log="output.log", output_cm=None, 
        report_time=10*u.picoseconds, sim_time=10*u.nanoseconds): 
    """
    Start and run an OpenMM NVT simulation with Langevin integrator at 2 fs 
    time step and 300 K. The cutoff distance for nonbonded interactions were 
    set at 1.2 nm and LJ switch distance at 1.0 nm, which commonly used with
    Charmm force field. Long-range nonbonded interactions were handled with PME.  

    Parameters
    ----------
    top_file : topology file (.top, .prmtop, ...)
        This is the topology file discribe all the interactions within the MD 
        system. 

    pdb_file : coordinates file (.gro, .pdb, ...)
        This is the molecule configuration file contains all the atom position
        and PBC (periodic boundary condition) box in the system. 
   
    check_point : None or check point file to load 
        
    GPU_index : Int or Str 
        The device # of GPU to use for running the simulation. Use Strings, '0,1'
        for example, to use more than 1 GPU
  
    output_traj : the trajectory file (.dcd)
        This is the file stores all the coordinates information of the MD 
        simulation results. 
  
    output_log : the log file (.log) 
        This file stores the MD simulation status, such as steps, time, potential
        energy, temperature, speed, etc.
 
    output_cm : the h5 file contains contact map information

    report_time : 10 ps
        The program writes its information to the output every 10 ps by default 

    sim_time : 10 ns
        The timespan of the simulation trajectory
    """
    
    top = pmd.load_file(top_file, xyz = pdb_file)
    
    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.2*u.nanometer,
                              switchDistance=1.0*u.nanometer, constraints=app.HBonds)
    dt = 0.002*u.picoseconds
    integrator = omm.LangevinIntegrator(300*u.kelvin, 1/u.picosecond, dt)
    
    try: 
        platform = omm.Platform_getPlatformByName("CUDA")
        properties = {'DeviceIndex': str(GPU_index), 'CudaPrecision': 'mixed'} 
    except Exception: 
        platform = omm.Platform_getPlatformByName("OpenCL")
        properties = {'DeviceIndex': str(GPU_index)} 
    
    simulation = app.Simulation(top.topology, system, integrator, platform, properties)
    
    simulation.context.setPositions(top.positions)
    
    simulation.minimizeEnergy()
    
    report_freq = int(report_time/dt)
    simulation.context.setVelocitiesToTemperature(10*u.kelvin, random.randint(1, 10000))
    simulation.reporters.append(app.DCDReporter(output_traj, report_freq))
    if output_cm: 
        simulation.reporters.append(ContactMapReporter(output_cm, report_freq))
    simulation.reporters.append(app.StateDataReporter(output_log,
            report_freq, step=True, time=True, speed=True,
            potentialEnergy=True, temperature=True, totalEnergy=True))
    simulation.reporters.append(app.CheckpointReporter('checkpnt.chk', report_freq))
    
    if check_point: 
        simulation.loadCheckpoint(check_point)
    nsteps = int(sim_time/dt)
    simulation.step(nsteps)
    
   
def openmm_simulate_amber_nvt(top_file, pdb_file, GPU_index=0, 
        output_traj="output.dcd", output_log="output.log", output_cm=None, 
        report_time=10*u.picoseconds, sim_time=10*u.nanoseconds): 
    """
    Start and run an OpenMM NVT simulation with Langevin integrator at 2 fs 
    time step and 300 K. The cutoff distance for nonbonded interactions were 
    set at 1.0 nm, which commonly used along with Amber force field. Long-range
    nonbonded interactions were handled with PME. 

    Parameters
    ----------
    top_file : topology file (.top, .prmtop, ...)
        This is the topology file discribe all the interactions within the MD 
        system. 

    pdb_file : coordinates file (.gro, .pdb, ...)
        This is the molecule configuration file contains all the atom position
        and PBC (periodic boundary condition) box in the system. 

    GPU_index : Int or Str 
        The device # of GPU to use for running the simulation. Use Strings, '0,1'
        for example, to use more than 1 GPU

    output_traj : the trajectory file (.dcd)
        This is the file stores all the coordinates information of the MD 
        simulation results. 

    output_log : the log file (.log) 
        This file stores the MD simulation status, such as steps, time, potential
        energy, temperature, speed, etc.

    report_time : 10 ps
        The program writes its information to the output every 10 ps by default 

    sim_time : 10 ns
        The timespan of the simulation trajectory
    """
    
    top = pmd.load_file(top_file, xyz = pdb_file)
    
    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.2*u.nanometer,
                              constraints=app.HBonds)
    dt = 0.002*u.picoseconds
    integrator = omm.LangevinIntegrator(300*u.kelvin, 1/u.picosecond, dt)
    
    try: 
        platform = omm.Platform_getPlatformByName("CUDA")
        properties = {'DeviceIndex': str(GPU_index), 'CudaPrecision': 'mixed'} 
    except Exception: 
        platform = omm.Platform_getPlatformByName("OpenCL")
        properties = {'DeviceIndex': str(GPU_index)} 
    
    simulation = app.Simulation(top.topology, system, integrator, platform, properties)
    
    simulation.context.setPositions(top.positions)
    
    simulation.minimizeEnergy()
    
    report_freq = int(report_time/dt)
    simulation.context.setVelocitiesToTemperature(10*u.kelvin, random.randint(1, 10000))
    simulation.reporters.append(app.DCDReporter(output_traj, report_freq))
    if output_cm:
        simulation.reporters.append(ContactMapReporter(output_cm, report_freq))
    simulation.reporters.append(app.StateDataReporter(output_log,
            report_freq, step=True, time=True, speed=True,
            potentialEnergy=True, temperature=True, totalEnergy=True))
    
    nsteps = int(sim_time/dt)
    simulation.step(nsteps)

#######################################################################################
def openmm_simulate_amber_fs_pep(pdb_file, top_file=None, check_point=None, GPU_index=0,
        output_traj="output.dcd", output_log="output.log", output_cm=None,
        report_time=10*u.picoseconds, sim_time=10*u.nanoseconds):
    """
    Start and run an OpenMM NVT simulation with Langevin integrator at 2 fs 
    time step and 300 K. The cutoff distance for nonbonded interactions were 
    set at 1.2 nm and LJ switch distance at 1.0 nm, which commonly used with
    Charmm force field. Long-range nonbonded interactions were handled with PME.  

    Parameters
    ----------
    pdb_file : coordinates file (.gro, .pdb, ...)
        This is the molecule configuration file contains all the atom position
        and PBC (periodic boundary condition) box in the system. 
   
    check_point : None or check point file to load 
        
    GPU_index : Int or Str 
        The device # of GPU to use for running the simulation. Use Strings, '0,1'
        for example, to use more than 1 GPU
  
    output_traj : the trajectory file (.dcd)
        This is the file stores all the coordinates information of the MD 
        simulation results. 
  
    output_log : the log file (.log) 
        This file stores the MD simulation status, such as steps, time, potential
        energy, temperature, speed, etc.
 
    output_cm : the h5 file contains contact map information

    report_time : 10 ps
        The program writes its information to the output every 10 ps by default 

    sim_time : 10 ns
        The timespan of the simulation trajectory
    """


    mytimer.mytime_label("pre_simulation",start=1)

    print("In openmm_simulation: 0");
    sys.stdout.flush()

    if top_file: 
        pdb = pmd.load_file(top_file, xyz = pdb_file)
        system = pdb.createSystem(nonbondedMethod=app.CutoffNonPeriodic, 
                nonbondedCutoff=1.0*u.nanometer, constraints=app.HBonds, 
                implicitSolvent=app.OBC1)
    else: 
        pdb = pmd.load_file(pdb_file)
        forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, 
                nonbondedCutoff=1.0*u.nanometer, constraints=app.HBonds)

    print("In openmm_simulation: 1");
    sys.stdout.flush()


    dt = 0.002*u.picoseconds
    integrator = omm.LangevinIntegrator(300*u.kelvin, 91.0/u.picosecond, dt)
    integrator.setConstraintTolerance(0.00001)

    print("In openmm_simulation: 2");
    sys.stdout.flush()


    try:
        platform = omm.Platform_getPlatformByName("CUDA")
        properties = {'DeviceIndex': str(GPU_index), 'CudaPrecision': 'mixed'}
    except Exception:
        platform = omm.Platform_getPlatformByName("OpenCL")
        properties = {'DeviceIndex': str(GPU_index)}

    print("In openmm_simulation: 3");
    sys.stdout.flush()


    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)

    print("In openmm_simulation: 4");
    sys.stdout.flush()


    simulation.context.setPositions(random.choice(pdb.get_coordinates())/10) #parmed \AA to OpenMM nm

    # equilibrate
    simulation.minimizeEnergy() 
    simulation.context.setVelocitiesToTemperature(300*u.kelvin, random.randint(1, 10000))
    simulation.step(int(100*u.picoseconds / (2*u.femtoseconds)))


    print("In openmm_simulation: 5");
    sys.stdout.flush()

    report_freq = int(report_time/dt)
    # simulation.reporters.append(app.DCDReporter(output_traj, report_freq))

    print("In openmm_simulation: 5a");
    sys.stdout.flush()

    if output_cm:
        simulation.reporters.append(ContactMapReporter(output_cm, report_freq))

    print("In openmm_simulation: 5b");
    sys.stdout.flush()

    '''
    simulation.reporters.append(app.StateDataReporter(output_log,
            report_freq, step=True, time=True, speed=True,
            potentialEnergy=True, temperature=True, totalEnergy=True))
    '''

    print("In openmm_simulation: 6");
    sys.stdout.flush()

    '''
    simulation.reporters.append(app.CheckpointReporter('checkpnt.chk', report_freq))
    '''

    print("In openmm_simulation: 6");
    sys.stdout.flush()

    '''
    if check_point:
        simulation.loadCheckpoint(check_point)
    '''
    nsteps = int(sim_time/dt)
    print(f'nsteps = {nsteps}, sim_time={sim_time}, dt = {dt}, report_time={report_time}')
    block = 0
    stop_file = 'stop.simulation'
    mytimer.mytime_label("pre_simulation",start=-1)

    '''
    while((not os.path.exists(stop_file)) and (not os.path.exists("../../../../aggregate/stop.aggregator"))):
          print(f'block={block} of the simulation')
          subprocess.getstatusoutput(f'touch start_block.{block}')
          sys.stdout.flush()
          mytimer.mytime_label("simulation.step",start=1)
          simulation.step(nsteps)
          mytop()
          mytimer.mytime_label("simulation.step",start=-1)
          subprocess.getstatusoutput(f'touch end_block.{block}')
    
          if(block == 5): #temporary, for debugging
              subprocess.getstatusoutput(f'touch {stop_file}')
    
          block += 1
    '''
    big_simulation = 0
    i3 = 0
    mytimer.mytime_label("simulation.big", start=1)              
    while(not os.path.exists("../../../../aggregate/stop.aggregator")):
        print(f'block={block} of the simulation')
        print(f'i3={i3}')
        # subprocess.getstatusoutput(f'touch start_block.{block}')
        sys.stdout.flush()
        if(not os.path.exists(stop_file)):
            mytimer.mytime_label("simulation.step",start=1)
            simulation.step(nsteps)
            mytop()
            mytimer.mytime_label("simulation.step",start=-1)
            # subprocess.getstatusoutput(f'touch end_block.{block}')
            block += 1
        else:
            mytimer.mytime_label("simulation.big", start=-1)
            mytimer.mytime_label("prepare_init", start=1)
            pdb_file = find_random_pdb(os.path.abspath("../../../../"))
            new_init = pmd.load_file(pdb_file).positions
            print(f"As initial condition using {pdb_file}")
            # print(type(new_init)); sys.stdout.flush()
            # print(new_init)
            # print(dir(new_init))
            simulation.context.setPositions(new_init)
            simulation.context.setVelocitiesToTemperature(300*u.kelvin, random.randint(1, 10000))
            subprocess.getstatusoutput(f"rm -f {stop_file}")
            mytimer.mytime_label("prepare_init", start=-1)
            mytimer.mytime_label("simulation.big", start=1)              
            big_simulation += 1
        if(i3 == 0):
            nsteps /= 5
        i3 += 1
    print("In openmm_simulation: 7");
    sys.stdout.flush()
#######################################################################################

def openmm_simulate_charmm_npt_z(top_file, pdb_file, check_point=None, GPU_index=0,
        output_traj="output.dcd", output_log="output.log", output_cm=None,
        report_time=10*u.picoseconds, sim_time=10*u.nanoseconds):
    """
    Start and run an OpenMM NVT simulation with Langevin integrator at 2 fs 
    time step and 300 K. The cutoff distance for nonbonded interactions were 
    set at 1.2 nm and LJ switch distance at 1.0 nm, which commonly used with
    Charmm force field. Long-range nonbonded interactions were handled with PME.  

    Parameters
    ----------
    top_file : topology file (.top, .prmtop, ...)
        This is the topology file discribe all the interactions within the MD 
        system. 

    pdb_file : coordinates file (.gro, .pdb, ...)
        This is the molecule configuration file contains all the atom position
        and PBC (periodic boundary condition) box in the system. 
   
    check_point : None or check point file to load 
        
    GPU_index : Int or Str 
        The device # of GPU to use for running the simulation. Use Strings, '0,1'
        for example, to use more than 1 GPU
  
    output_traj : the trajectory file (.dcd)
        This is the file stores all the coordinates information of the MD 
        simulation results. 
  
    output_log : the log file (.log) 
        This file stores the MD simulation status, such as steps, time, potential
        energy, temperature, speed, etc.
 
    output_cm : the h5 file contains contact map information

    report_time : 10 ps
        The program writes its information to the output every 10 ps by default 

    sim_time : 10 ns
        The timespan of the simulation trajectory
    """

    top = pmd.load_file(top_file, xyz = pdb_file)

    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.2*u.nanometer,
                              switchDistance=1.0*u.nanometer, constraints=app.HBonds)
    dt = 0.002*u.picoseconds
    integrator = omm.LangevinIntegrator(300*u.kelvin, 1/u.picosecond, dt)

    system.addForce(omm.MonteCarloAnisotropicBarostat((1, 1, 1)*u.bar, 300*u.kelvin, False, False, True)) 

    try:
        platform = omm.Platform_getPlatformByName("CUDA")
        properties = {'DeviceIndex': str(GPU_index), 'CudaPrecision': 'mixed'}
    except Exception:
        platform = omm.Platform_getPlatformByName("OpenCL")
        properties = {'DeviceIndex': str(GPU_index)}

    simulation = app.Simulation(top.topology, system, integrator, platform, properties)

    simulation.context.setPositions(top.positions)

    simulation.minimizeEnergy()

    report_freq = int(report_time/dt)
    simulation.context.setVelocitiesToTemperature(10*u.kelvin, random.randint(1, 10000))
    simulation.reporters.append(app.DCDReporter(output_traj, report_freq))
    if output_cm:
        simulation.reporters.append(ContactMapReporter(output_cm, report_freq))
    simulation.reporters.append(app.StateDataReporter(output_log,
            report_freq, step=True, time=True, speed=True,
            potentialEnergy=True, temperature=True, totalEnergy=True))
    simulation.reporters.append(app.CheckpointReporter('checkpnt.chk', report_freq))

    if check_point:
        simulation.loadCheckpoint(check_point)
    nsteps = int(sim_time/dt)
    simulation.step(nsteps)


def openmm_simulate_amber_npt(top_file, pdb_file, check_point, GPU_index=0,
        output_traj="output.dcd", output_log="output.log", output_cm=None,
        report_time=10*u.picoseconds, sim_time=10*u.nanoseconds):
    """
    Start and run an OpenMM NVT simulation with Langevin integrator at 2 fs 
    time step and 300 K. The cutoff distance for nonbonded interactions were 
    set at 1.0 nm, which commonly used along with Amber force field. Long-range
    nonbonded interactions were handled with PME. 

    Parameters
    ----------
    top_file : topology file (.top, .prmtop, ...)
        This is the topology file discribe all the interactions within the MD 
        system. 

    pdb_file : coordinates file (.gro, .pdb, ...)
        This is the molecule configuration file contains all the atom position
        and PBC (periodic boundary condition) box in the system. 

    GPU_index : Int or Str 
        The device # of GPU to use for running the simulation. Use Strings, '0,1'
        for example, to use more than 1 GPU

    output_traj : the trajectory file (.dcd)
        This is the file stores all the coordinates information of the MD 
        simulation results. 

    output_log : the log file (.log) 
        This file stores the MD simulation status, such as steps, time, potential
        energy, temperature, speed, etc.

    report_time : 10 ps
        The program writes its information to the output every 10 ps by default 

    sim_time : 10 ns
        The timespan of the simulation trajectory
    """

    top = pmd.load_file(top_file, xyz = pdb_file)

    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*u.nanometer,
                              constraints=app.HBonds)
    dt = 0.002*u.picoseconds
    integrator = omm.LangevinIntegrator(300*u.kelvin, 1/u.picosecond, dt)
    system.addForce(omm.MonteCarloBarostat(1*u.bar, 300*u.kelvin))

    try:
        platform = omm.Platform_getPlatformByName("CUDA")
        properties = {'DeviceIndex': str(GPU_index), 'CudaPrecision': 'mixed'}
    except Exception:
        platform = omm.Platform_getPlatformByName("OpenCL")
        properties = {'DeviceIndex': str(GPU_index)}

    simulation = app.Simulation(top.topology, system, integrator, platform, properties)

    simulation.context.setPositions(top.positions)

    simulation.minimizeEnergy()

    report_freq = int(report_time/dt)
    simulation.context.setVelocitiesToTemperature(10*u.kelvin, random.randint(1, 10000))
    simulation.reporters.append(app.DCDReporter(output_traj, report_freq))
    if output_cm:
        simulation.reporters.append(ContactMapReporter(output_cm, report_freq))
    simulation.reporters.append(app.StateDataReporter(output_log,
            report_freq, step=True, time=True, speed=True,
            potentialEnergy=True, temperature=True, totalEnergy=True))
    simulation.reporters.append(app.CheckpointReporter('checkpnt.chk', report_freq))

    if check_point:
        simulation.loadCheckpoint(check_point)
    nsteps = int(sim_time/dt)
    simulation.step(nsteps)
