## Downloading data generated from this study:

mkdir all_counts/
cd all_counts/

wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6819/E-MTAB-6819.sdrf.txt --no-check-certificate

wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6819/E-MTAB-6819.processed.1.zip --no-check-certificate
unzip E-MTAB-6819.processed.1.zip
rm E-MTAB-6819.processed.1.zip

cd -

## Downloading data from other studies for remapping:

mkdir other_counts/
cd other_counts/

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100597/suppl/GSE100597%5Fcount%5Ftable%5FQC%5Ffiltered%2Etxt%2Egz # Mohammed

wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3929/E-MTAB-3929.processed.1.zip --no-check-certificate # Petropoulos
mkdir E-MTAB-3929
unzip E-MTAB-3929.processed.1.zip -d E-MTAB-3929
rm E-MTAB-3929.processed.1.zip

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74767/suppl/GSE74767%5FSC3seq%5FCy%5FProcessedData%2Etxt%2Egz # Nakamura

cd -

## Downloading the bulk RNA-seq data

mkdir bulk_counts/
cd bulk_counts/

wget https://jmlab-gitlab.cruk.cam.ac.uk/aaron/NaiveHESC2016-DataFiles/raw/master/bulk/Guo2017_rawReadcountRNAseq.txt
wget https://jmlab-gitlab.cruk.cam.ac.uk/aaron/NaiveHESC2016-DataFiles/raw/master/bulk/Pastor2016_rawReadcountRNAseq.txt.gz
wget https://jmlab-gitlab.cruk.cam.ac.uk/aaron/NaiveHESC2016-DataFiles/raw/master/bulk/Theunissen2016_rawReadcountRNAseq.txt.gz

cd -

## Also cloning the tools (using the last commit known to work).
#git clone https://github.com/LTLA/CRUKTools tools
#cd tools
#git reset --hard 915706ac016f6f59e74b4ec266f06349293b997f
