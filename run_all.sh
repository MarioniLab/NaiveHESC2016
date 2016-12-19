set -e
set -u
echo "knitr::knit('preprocess.Rmd')" | R --no-save --vanilla
echo "knitr::knit('overall_analysis.Rmd')" | R --no-save --vanilla
echo "knitr::knit('naive_analysis.Rmd')" | R --no-save --vanilla
echo "knitr::knit('primed_analysis.Rmd')" | R --no-save --vanilla
echo "knitr::knit('remapping.Rmd')" | R --no-save --vanilla

