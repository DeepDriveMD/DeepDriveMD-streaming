ip='129.114.17.185'
pass="KpdsbAJD9zDH"
#old='129.114.17.233'
export RMQ_HOSTNAME=$ip
export RMQ_PORT=5672
export RMQ_USERNAME="iyakushin"
export RMQ_PASSWORD=$pass
export RADICAL_PILOT_DBURL=mongodb://iyakushin:${pass}@${ip}:27017/entk-test
export RADICAL_PROFILE=True
export RADICAL_PILOT_PROFILE=True
export RADICAL_ENTK_PROFILE=True

export PYTHONPATH=/gpfs/alpine/scratch/iyakushin/csc299/Test3/entk_cvae_md_devel3/entk_cvae_md/Outlier_search:$PYTHONPATH

python adios2pandasR5.py


