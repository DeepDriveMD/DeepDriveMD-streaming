dir=/p/gpfs1/yakushin/radical.pilot.sandbox/$(ls -d ../re.* | cut -d"/" -f2)/pilot.0000
echo $dir
grep intersect ${dir}/unit.00*/*.out
