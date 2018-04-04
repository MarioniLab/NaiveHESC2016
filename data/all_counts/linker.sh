# Aggregates all of the count files into a single directory, appropriately renamed.

for x in $(ls .. | grep "real-")
do 
    run=$(basename $x | sed -E "s/real-([0-9]+)_.*/\1/")
    ln ../$x/genic_counts.tsv genic_counts_$run.tsv
done
