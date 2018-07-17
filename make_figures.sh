Rscript Figure1.R
Rscript Figure2.R
Rscript Figure3.R
Rscript Figure4.R



targetdir=~/Dropbox/primed_vs_naive_ESC_analysis/latest_analysis/
if [ -e $targetdir ]
then
    cp *.html $targetdir
    cp -r results-* $targetdir
    cp -r Figures $targetdir
else 
    echo "'$targetdir' does not exist!"
    exit 1    
fi

