# assigning peptides to taxonomic groups and outputting formatted csvs
# only for the data that uses TFG-based database, because the annotation
# files are in excel and functions used are unique to this.

Rscript peptide-assignment-post-mcl.R "../data/tfg-all-database/" "tfg-all" &
Rscript peptide-assignment-post-mcl.R "../data/tfg-t0-database/" "tfg-t0" &

