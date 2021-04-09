python timers.py

./120run.sh
./120test.sh > 120test.csv


#python outliers.py

cat 120test.csv

rm -f *.gz
tar zcvf out.tar.gz *.csv
