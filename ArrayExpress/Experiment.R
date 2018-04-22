all.out <- list()
all.out[["MAGE-TAB Version"]] <- "1.1"
all.out[["Investigation Title"]] <- "Single-cell RNA-seq of naive and primed human embryonic stem cells"
all.out[["Experiment Description"]] <- "This study aims to profile the transcriptomes of single naive and primed human embryonic stem cells. Cells from the H9 line were cultured to select for naive or primed phenotypes, and a sequencing library was generated from each single cell using the Smart-seq2 method. This was repeated for multiple experimental batches, i.e., independent cultures. Batch 1 consists of sequencing runs 2383 and 2384; batch 2 consists of runs 2678, 2679, 2739 and 2740; and batch 3 consists of runs 2780 and 2781. Transcriptional profiles for all cells were obtained from the sequencing data and used to explore substructure and heterogeneity in the population for each phenotype."

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

all.out[["Protocol Description"]] <- c("Human H9-NK2 ESCs were kindly provided by Austin Smith." ,
                                       "Naive hESCs were grown in 6-well dishes on mouse embryonic fibroblasts in N2B27 supplemented with human LIF, 1 uM Chiron, 1 uM PD03 and 2 uM Go6983. One passage before sorting, cells were plated on 6-well plates coated with Matrigel (growth-factor reduced). Primed hESCs were grown in 6-well dishes coated with Vitronectin in E8 media.",
                                       "hESCs were dissociated with Accutase and sorted with a BD Aria Cell sorter, gating for cell size and granularity. Single-cells were sorted in 2uL of Lysis Buffer (0.2% v/v Triton X-100 (Sigma-Aldrich, cat. no. T9284) with 2U/ul RNase Inhibitor (Clontech, cat. no. 2313A)) in 96 well plates, spun down and immediately frozen at -80 degrees Celsius.",
                                       "The cDNA libraries for sequencing were prepared using Nextera XT DNA Sample Preparation Kit (Illumina, cat. no. FC-131-1096), according to the protocol supplied by Fluidigm (PN 100-5950 B1). Libraries from 96 single cells were pooled and purified using AMPure XP beads (Beckman Coulter).",
                                       "Pooled samples were sequenced on an Illumina HiSeq 2500 instrument, using paired-end 100-bp reads.",
                                       "Reads were aligned to the hg38 build of the human genome (with additional ERCC sequences) using subread v1.6.1. in paired-end RNA-seq mode with unique mapping. The number of read pairs mapped to the exonic regions of each gene was then counted for each library, using the featureCounts function in Rsubread v1.28.1 with Ensembl GRCh38 version 91. Only alignments with mapping quality scores above 10 were considered during counting.",
                                        "A cDNA library was prepared from each sorted single cell following the SmartSeq2 protocol. Briefly, oligo-dT primer, dNTPs (ThermoFisher, cat. no. 10319879) and ERCC RNA Spike-In Mix (1:25,000,000 final dilution, Ambion, cat. no. 4456740) were added to the single-cell lysates, and reverse transcription and PCR were performed.")
all.out[["Protocol Hardware"]] <- c("",
                                    "",
                                    "BD FACSAria cell sorter", 
                                    "",
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
all.out[["Comment[AdditionalFile:fa]"]] <- "spikes.fa"
all.out[["Comment[AdditionalFile:txt]"]] <- "spikes.txt"

unlink("idf.tsv")
for (x in names(all.out)) {
    write(file="idf.tsv", paste0(c(x, all.out[[x]]), collapse="\t"), append=TRUE)
}

##########################################################################3
# Computing the spike-in quantity per well.

# ERCC data taken from https://www.thermofisher.com/order/catalog/product/4456740,
# under the link "ERCC Controls Analysis: ERCC RNA Spike-In Control Mixes (English)".
ercc.dil <- 25e6
ercc.data <- read.table("cms_095046.txt", header=TRUE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)
ercc.id <- ercc.data[,"ERCC ID"]
ercc.quant <- ercc.data[,"concentration in Mix 1 (attomoles/ul)"] / ercc.dil

spike.dir <- "spike-data"
dir.create(spike.dir, showWarnings=FALSE)
write.table(file=file.path(spike.dir, "spikes.txt"), sep="\t", quote=FALSE, row.names=FALSE,
    data.frame(Name=ercc.id, "Attomole/well"=ercc.quant, check.names=FALSE))


