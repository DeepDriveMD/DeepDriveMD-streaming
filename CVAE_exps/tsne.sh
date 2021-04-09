lastN=$1
python best.py
python project.py $lastN
python rmsd.py 40 $lastN
