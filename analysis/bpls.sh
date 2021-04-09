for bpf in $(ls -d ../aggregate/agg*.bp); do
    echo $bpf
    bpls $bpf
    echo "=========="
done
