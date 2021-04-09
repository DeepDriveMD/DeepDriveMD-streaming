import glob
import subprocess
import pandas as pd
import os

session = subprocess.getstatusoutput('ls -d ../re.* | cut -d "/" -f2')[1]
dir = f'/p/gpfs1/yakushin/radical.pilot.sandbox/{session}/pilot.0000'

stdouts = glob.glob(f'{dir}/unit*/*.out')

pf = pd.DataFrame(columns = ["label", "start", "gps", "date", "file", "line", "unit", "time"])

labels = []
starts = []
gpss = []
dates = []
files = []
nlines = []
units = []
times = []



fn = f"{session}.csv"

for s in stdouts:
    unit = int(s.split("/")[-2].replace("unit.",""))
    with open(s) as f:
        lines = list(filter(lambda x: x.find("TLaBeL") == 0 and x.find("Testing") == -1, f.readlines()))
        for line in lines:
            tokens = line.split("|")
            labels.append(tokens[1])
            starts.append(int(tokens[2]))
            gpss.append(int(float(tokens[3])))
            dates.append(tokens[4])
            files.append(os.path.basename(tokens[5]))
            nlines.append(int(tokens[6].strip()))
            units.append(unit)
            times.append(float(tokens[7].strip()))


print(len(labels))


pf['label'] = labels
pf['start'] = starts
pf['gps'] = gpss
pf['date'] = dates
pf['file'] = files
pf['line'] = nlines
pf['unit'] = units
pf['time'] = times

pf.to_csv(fn)


'''
print(subprocess.getstatusoutput('/gpfs/alpine/scratch/iyakushin/csc299/radical.pilot.sandbox/$(ls -tr /gpfs/alpine/scratch/iyakushin/csc299/radical.pilot.sandbox/ | tail -1)/pilot.0000'))

dir=/gpfs/alpine/scratch/iyakushin/csc299/radical.pilot.sandbox/$(ls -tr /gpfs/alpine/scratch/iyakushin/csc299/radical.pilot.sandbox/ | tail -1)/pilot.0000
echo $dir
grep TLaBeL ${dir}/unit.0000*/STDOUT | grep -v Testing

['TLaBeL', 'openmm_simulate', '1', '1599750189.0', 'Thu Sep 10 11:03:09 2020', '/gpfs/alpine/csc299/scratch/iyakushin/Test3/entk_cvae_md/MD_exps/fs-pep/run_openmm.py', '104\n']




ignore_index = True

'''
