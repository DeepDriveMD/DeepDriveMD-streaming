dir=/p/gpfs1/yakushin/radical.pilot.sandbox/$(ls -d ../re* | cut -d"/" -f2)/pilot.0000
echo $dir
grep initial ${dir}/unit.*/*.out | perl -ne '@a=split(/\//);print $a[-1]' | sort | uniq -c
