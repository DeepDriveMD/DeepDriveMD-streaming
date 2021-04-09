dir=/gpfs/alpine/scratch/iyakushin/csc299/radical.pilot.sandbox/$(ls -d ../re.* | cut -d'/' -f2)/pilot.0000
echo $dir
grep TLaBeL ${dir}/unit.0000*/STDOUT | grep -v Testing
