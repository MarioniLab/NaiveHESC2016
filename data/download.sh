## Downloading data generated from this study:

cd all_counts/
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6819/E-MTAB-6819.sdrf.txt
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6819/E-MTAB-6819.processed.1.zip
unzip E-MTAB-6819.processed.1.zip
rm E-MTAB-6819.processed.1.zip
cd -

## Downloading data from other studies for remapping:

cd other_counts/

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100597/suppl/GSE100597%5Fcount%5Ftable%5FQC%5Ffiltered%2Etxt%2Egz # Mohammed

wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3929/E-MTAB-3929.processed.1.zip # Petropoulos
mkdir E-MTAB-3929
unzip E-MTAB-3929.processed.1.zip -d E-MTAB-3929
rm E-MTAB-3929.processed.1.zip

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74767/suppl/GSE74767%5FSC3seq%5FCy%5FProcessedData%2Etxt%2Egz # Nakamura

cd -

## Also cloning the tools (using the last commit known to work).
#git clone https://github.com/LTLA/CRUKTools tools
#cd tools
#git reset --hard 915706ac016f6f59e74b4ec266f06349293b997f
