dir=../MD_exps/fs-pep/
running=$dir/running
stopped=$dir/stopped
new=$dir/new
all=$dir/all


echo "Running"
ls -1 $running | wc

ls -1 $running
echo "============="

echo "Stopped"
ls -1 $stopped | wc

echo "============="

echo "New"
ls -1 $new | wc

ls -1 $new
echo "============="

echo "All"
ls -1 $all | wc
