set -e
set -u

echo "rmarkdown::render('analysis/preprocess.Rmd', output_file = 'preprocess.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/overall_analysis.Rmd', output_file = 'overall_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/naive_analysis.Rmd', output_file = 'naive_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/transition_analysis.Rmd', output_file = 'transition_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/primed_analysis.Rmd', output_file = 'primed_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/diff.Rmd', output_file = 'diff.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/correlations.Rmd', output_file = 'correlations.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/remapping.Rmd', output_file = 'remapping.html')" | R --no-save --vanilla
echo "rmarkdown::render('analysis/qualitycheck.Rmd', output_file = 'qualitycheck.html')" | R --no-save --vanilla

targetdir="~/Desktop/Ergebnisse
"
if [ -e $targetdir ]
then
    cp *.html $targetdir
    cp -r results-* $targetdir
    cp -r Figures $targetdir
else 
    echo "'$targetdir' does not exist!"
    exit 1    
fi