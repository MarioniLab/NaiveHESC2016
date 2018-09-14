set -e
set -u

echo "rmarkdown::render('preprocess.Rmd', output_file = 'preprocess.html')" | R --no-save --vanilla
echo "rmarkdown::render('overall_analysis.Rmd', output_file = 'overall_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('naive_analysis.Rmd', output_file = 'naive_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('transition_analysis.Rmd', output_file = 'transition_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('primed_analysis.Rmd', output_file = 'primed_analysis.html')" | R --no-save --vanilla
echo "rmarkdown::render('diff.Rmd', output_file = 'diff.html')" | R --no-save --vanilla
echo "rmarkdown::render('correlations.Rmd', output_file = 'correlations.html')" | R --no-save --vanilla
echo "rmarkdown::render('remapping.Rmd', output_file = 'remapping.html')" | R --no-save --vanilla
echo "rmarkdown::render('qualitycheck.Rmd', output_file = 'qualitycheck.html')" | R --no-save --vanilla

