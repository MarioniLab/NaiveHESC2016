all.out <- list()
all.out[["MAGE-TAB Version"]] <- "1.1"
all.out[["Investigation Title"]] <- "Single-cell RNA-seq of naive and primed human embryonic stem cells"
all.out[["Experiment Description"]] <- "This study aims to profile the transcriptomes of single naive and primed human embryonic stem cells. Cells from the H9 line were cultured to select for naive or primed phenotypes, and a sequencing library was generated from each single cell using the Smart-seq2 method. This was performed across multiple experimental batches (i.e., independent cultures) and high-throughput sequencing runs. Transcriptional profiles for all cells were obtained from the sequencing data and used to explore substructure and heterogeneity in the population for each phenotype."

all.out[["Experimental Design"]] <- c("cell type comparison design")
all.out[["Experimental Design Term Source REF"]] <- c("EFO") 
all.out[["Experimental Design Term Accession Number"]] <- c("EFO_0001745")

all.out[["Experimental Factor Name"]] <- "phenotype"
all.out[["Experimental Factor Type"]] <- "phenotype"
all.out[["Experimental Factor Term Source REF"]] <- ""
all.out[["Experimental Factor Term Accession Number"]] <- ""

all.out[["Person Last Name"]] <- "Lun"
all.out[["Person First Name"]] <- "Aaron"
all.out[["Person Mid Initials"]] <- "TL"
all.out[["Person Email"]] <- "aaron.lun@cruk.cam.ac.uk"
all.out[["Person Phone"]] <- ""
all.out[["Person Fax"]] <- ""      
all.out[["Person Address"]] <- "University of Cambridge Li Ka Shing Centre Robinson Way Cambridge CB2 0RE United Kingdom"
all.out[["Person Affiliation"]] <- "Cancer Research UK Cambridge Institute"
all.out[["Person Roles"]] <- "submitter"

all.out[["Protocol Name"]] <- c("Obtaining H9 cells",
                                "Culturing H9 cells",
                                "Extracting RNA",
                                "Creating libraries",
                                "Sequencing libraries",
                                "Assigning reads to genes",
                                "Reverse transcription")
all.out[["Protocol Type"]] <- c("sample collection protocol",
                                "growth protocol",
                                "nucleic acid extraction protocol",
                                "nucleic acid library construction protocol",
                                "nucleic acid sequencing protocol",
                                "high throughput sequence alignment protocol",
                                "conversion protocol")
all.out[["Protocol Term Source REF"]] <- c("EFO", 
                                           "EFO",   
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO")
all.out[["Protocol Term Accession Number"]] <- c("EFO_0005518",
                                                 "EFO_0003789",
                                                 "EFO_0002944",
                                                 "EFO_0004184",
                                                 "EFO_0004170",
                                                 "EFO_0004917",
                                                 "EFO_0005520")

all.out[["Protocol Description"]] <- c("HeLa cells were obtained from American Type Culture Collection. To generate clones expressing dCas9-KRAB, HeLa cells were transduced with lentivirus containing the pHR-SFFV-dCAS9-BFP-KRAB vector along with polybrene (5 ug/ml, Sigma), incubated for 72 hours with medium replacement after 24 hours, and sorted for the BFP-expressing cells using a BD FACSAria III cell sorter.",
                                       "HeLa cells were maintained in Dulbecco's modified Eagle's medium (Sigma Aldrich, D6429) supplemented with 10% fetal bovine serum (Thermo Fisher Scientific) and cultured at 37 degrees Celsius with with 5% CO2.",
                                       "RNA (1 ug) was extracted with the RNeasy Kit (QIAGEN) and treated with DNase I following the manufacturer's instructions. RNA quality was assessed using a Total RNA Nano chip with a 2100 Bioanalyzer instrument (Agilent).",
                                       "RNA-seq libraries were prepared from HeLA cells using TruSeq Stranded Total RNA Kit with Ribo-Zero Gold (Illumina, RS-122-2303). Library quality was assessed using a DNA1000 chip with a 2100 Bioanalyzer instrument (Agilent).",
                                       "Indexed libraries were PCR amplified and sequenced on multiple lanes of an Illumina Hiseq 2500 instrument to obtain 125 bp paired-end reads.",
                                       "Reads were aligned to the hg38 build of the human genome using subread v1.5.0 in paired-end RNA-seq mode with unique mapping. The number of read pairs mapped to the exonic regions of each gene was then counted for each library, using the featureCounts function in Rsubread v1.22.3 with Ensembl GRCh38 version 83. Only alignments with mapping quality scores above 10 were considered during counting.",
                                       "The QuantiTect Reverse Transcription Kit (QIAGEN) was used for cDNA synthesis including an additional step to eliminate genomic DNA contamination.")
all.out[["Protocol Hardware"]] <- c("BD FACSAria III cell sorter",
                                    "", 
                                    "2100 Bioanalyzer",
                                    "2100 Bioanalyzer",
                                    "Illumina Hiseq 2500",
                                    "",
                                    "")
all.out[["Protocol Software"]] <- c("", 
                                    "",
                                    "",
                                    "",
                                    "",
                                    "(R)subread",
                                    "")
                                 
all.out[["Term Source Name"]] <- "EFO"
all.out[["Term Source File"]] <- "http://www.ebi.ac.uk/efo/"
all.out[["Term Source Version"]] <- ""
all.out[["SDRF File"]] <- "sdrf.tsv"
all.out[["Public Release Date"]] <- "2017-05-31"
all.out[["Comment[AEExperimentType]"]] <- "RNA-seq of coding RNA"

unlink("idf.tsv")
for (x in names(all.out)) {
    write(file="idf.tsv", paste0(c(x, all.out[[x]]), collapse="\t"), append=TRUE)
}


