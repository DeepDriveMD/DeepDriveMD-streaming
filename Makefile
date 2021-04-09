test7:
	./run_test7.sh
test6:
	./run_test6.sh
test5:
	./run_test5.sh
test4:
	./run_test4.sh
test3a:
	./run_test3a.sh
test3:
	./run_test3.sh
test2:
	./run_test2.sh
test:
	./run_test.sh
run:
	./run1.sh
clean:
	rm -rf MD_exps/fs-pep/{all,running,stopped,new} MD_to_CVAE/cvae_input.* CVAE_exps/cvae_runs_*  *.prof *~ __pycache__ */__pycache__ */*/__pycache__ */*~ */*/*~ MD_to_CVAE/cvae_input.bp Outlier_search/{outlier_pdbs,tmp} Outlier_search/{*.json,*.lock,*.pickle} aggregate/{stop.aggregator,aggregator.log,agg*.bp,aggregate.cprofile}  analysis/*.{out,err,csv,gz,tar} Outlier_search/archive *.log re.* ./Outlier_search/__pycache__  ./MD_exps/MD_utils/__pycache__  ./CVAE_exps/cvae/__pycache__ ./aggregate/__pycache__

clean_test:
	rm -rf  re.session.login*.iyakushin* *.prof *~ __pycache__ */__pycache__ */*/__pycache__ */*~ */*/*~ MD_to_CVAE/cvae_input.bp Outlier_search/{outlier_pdbs,tmp,*.lock,*.pickle} Outlier_search/*.json MD_exps/fs-pep/{all,new,running,stopped} aggregate/{stop.aggregator,aggregator.log,aggregator.bp}
