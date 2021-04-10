# DeepDriveMD-S: Streaming version of DeepDriveMD

* To set up the environment for Radical ENTK to run the program:
  ```
  source radical.sh
  ```
* To run:
  ```
  make
  ```	
* One might need to change various paths in `summit_md_test7.py` to run on a different account.
* The environment for each Task is set by sourcing `powerai.sh` that sets up a python environment and loads ADIOS2 module.
* This python environment is build on top of powerai python environment from IBM. Various MD packages like OpenMM, parmed are added.
* The program uses TensorFlow, ADIOS2, OpenMM, parmed, RAPIDS' GPU-accelearted implementation of DBSCAN, LockFile besides standard python packages.
