export RMQ_HOSTNAME=cz-deepdrivemd.apps.czapps.llnl.gov
export RMQ_PASSWORD=x4ESqyNyQ1F636qjVfLmmpANWiBra/V2NMfiBEDd5IklzzF7IlONbwMQlPWgpqz6
export RMQ_PORT=31448
export RMQ_USERNAME=deepdrive
export RMQ_VHOST=/deepdrive
mongodb_host=cz-deepdrivemd.apps.czapps.llnl.gov
mongodb_port=31945
mongodb_database=deepdriveDB
mongodb_user=deepdrive
mongodb_password=G9eHTINgrFOU6IgP0gPzCMDi9vlGslows91qzbXwvPOgqmv674VjDTCjn3%2FxXAo2
export RADICAL_PILOT_DBURL="mongodb+ssl://$mongodb_user:$mongodb_password@$mongodb_host:$mongodb_port/$mongodb_database?tlsAllowInvalidCertificates=True"


export RADICAL_PROFILE=True
export RADICAL_PILOT_PROFILE=True
export RADICAL_ENTK_PROFILE=True
export RADICAL_LOG_LVL=DEBUG
export RADICAL_LOG_TGT=radical.log
export RMQ_SSL=True


#export PYTHONPATH=$GPFS1/DDMD/26/entk_cvae_md/Outlier_search:$PYTHONPATH
#export PYTHON=`which python`

#export RMQ_HOSTNAME=129.114.17.185
#export RMQ_PORT=5672
#export RMQ_USERNAME=hyperrct
#export RMQ_PASSWORD=h1p3rrc7
#unset RMQ_SSL
#unset RMQ_VHOST



python summit_md_test7.py

