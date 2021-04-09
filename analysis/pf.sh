echo "=============="
hostname
source /ccs/home/iyakushin/etc/openmm5.sh
which python
cd /gpfs/alpine/scratch/iyakushin/csc299/Test3/entk_cvae_md_devel3/entk_cvae_md/analysis
python adios2pandasR4.py $1 $2
echo "=============="
