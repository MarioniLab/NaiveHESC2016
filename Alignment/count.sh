set -u  

for dir in real-2383_20161024 real-2384_20161024 real-2677_20161024 real-2678_20161024 real-2739_20161024 real-2740_20161024 real-2780_20161024 real-2781_20161024

do
	cd /lustre/jmlab/messme01/vMeyenn/$dir
	bsub -J "readcount_1" -R "rusage[mem=16000]" -n 1 -o std.out -e std.err "Rscript ./counter.R" 
	cd -	
done