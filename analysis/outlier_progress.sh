for outlier in $(ls outlier*.csv | sort -t '_' -k 2 -g); do
    echo $outlier
    best=`cut -d',' -f3 $outlier | grep -v R | sort -gr | tail -1`
    echo $best
done
