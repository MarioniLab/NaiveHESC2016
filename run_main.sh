set -e
set -u
echo "knitr::knit('preprocess.Rmd')" | Rdevel --no-save --vanilla
echo "knitr::knit('overall_analysis.Rmd')" | Rdevel --no-save --vanilla
echo "knitr::knit('naive_analysis.Rmd')" | Rdevel --no-save --vanilla
echo "knitr::knit('transition_analysis.Rmd')" | Rdevel --no-save --vanilla
echo "knitr::knit('primed_analysis.Rmd')" | Rdevel --no-save --vanilla
echo "knitr::knit('diff.Rmd')" | Rdevel --no-save --vanilla
echo "knitr::knit('correlations.Rmd')" | Rdevel --no-save --vanilla
echo "knitr::knit('remapping.Rmd')" | Rdevel --no-save --vanilla

echo "rmarkdown::render('preprocess.md')" | R --no-save --vanilla
echo "rmarkdown::render('overall_analysis.md')" | R --no-save --vanilla
echo "rmarkdown::render('naive_analysis.md')" | R --no-save --vanilla
echo "rmarkdown::render('transition_analysis.md')" | R --no-save --vanilla
echo "rmarkdown::render('primed_analysis.md')" | R --no-save --vanilla
echo "rmarkdown::render('diff.md')" | R --no-save --vanilla
echo "rmarkdown::render('correlations.md')" | R --no-save --vanilla
echo "rmarkdown::render('remapping.md')" | R --no-save --vanilla

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

