cd data/counts/
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6819/E-MTAB-6819.sdrf.txt
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6819/E-MTAB-6819.processed.1.zip
unzip E-MTAB-6819.processed.1.zip
rm E-MTAB-6819.processed.1.zip
cd -


# Also cloning the tools (using the last commit known to work).
#git clone https://github.com/LTLA/CRUKTools tools
#cd tools
#git reset --hard 915706ac016f6f59e74b4ec266f06349293b997f