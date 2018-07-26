# This script generates new Rmarkdown files using only G1-related cells 
# in naive_analysis.Rmd and primed_analysis.Rmd. The idea is to check
# that the same conclusions hold when the effect of cell cycle is removed.

for mode in naive primed
do
    newfile=analysis/g1_${mode}_analysis.Rmd
    cat ${mode}_analysis.Rmd | sed "s/figure-${mode}/figure-g1_${mode}/g" \
        | sed -r 's/(sce <- readRDS\("sce_all.rds"\))/\1\nsce <- sce[,sce$phase=="G1"]/' \
        | sed "s/resdir <- \"results-${mode}\"/resdir <- \"results-g1_${mode}\"/" \
        | sed "s/saveRDS(sce_${mode}, file=\"sce_${mode}.rds\")//" \
        | sed "s/embryonic stem cells: ${mode}-only analysis/embryonic stem cells: G1-only, ${mode}-only analysis/" \
        | sed -r "s/(This is more sensitive than doing so on the entire data set, where the naive\/primed differences dominate.)/\1\nWe subset to only use G1 cells, to check that the cell cycle phase does not have any effect on the analysis./" \
        | sed "s/We save our object for later use and/We/" > ${newfile}

    if [[ ${mode} == "naive" ]]
    then
        top="We also check the cell cycle phase for each of the three clusters."
        bottom="We also examine the location of cluster 3 on the overall PCA plot."
        sed "/^${top}/,/^${bottom}/ { /^${bottom}/b; d }" ${newfile} > .tmp.Rmd
    else 
        top="### Effect of the cell cycle"
        bottom="# Wrapping up"
        sed "/^${top}/,/^${bottom}/ { /^${bottom}/b; d }" ${newfile} > .tmp.Rmd
    fi
    mv .tmp.Rmd ${newfile}
 done



    
    
