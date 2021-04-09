export PYTHONPATH=../Outlier_search/:../misc/:$PYTHONPATH

python best.py
python project.py 10000
python rmsd.py 40 10000
python embeddings.py 40 3
python embeddings.py 40 2
