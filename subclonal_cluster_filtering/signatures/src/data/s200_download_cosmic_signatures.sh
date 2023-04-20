#Downloads COSMIC v2 signatures and fixes the header and last line (missing
#\n).
#Also changes the format to that of SigProfiler signature tables, i.e. with
# - 'Type'      mutating pyrimidine
# - 'SubType'   trinucleotide context
#This means the 'Somatic_Mutation_Type' (e.g. 'A[C>A]A') column is dropped

#File name variables
    IN_FILE="http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
    OUT_DIR="data/raw/reference_signatures/"
    OUT_FILE="${OUT_DIR}COSMIC_signatures.tsv"

#Make sure the output directory exists
    mkdir -p $OUT_DIR

#Download
    wget -nv $IN_FILE -O ${OUT_FILE}.temp

#Fix header and rename file
    head -n 1 ${OUT_FILE}.temp |
        cut -f4- |
        perl -ne '$_ =~ s/ /_/g; print "Type\tSubType\t$_"' \
        > $OUT_FILE

#Append the rest of the table and fix the last line
    tail -n +2 ${OUT_FILE}.temp |
        cut -f-2,4- |
        perl -ne 'chomp; print "$_\n";' |
        sort -k 1 \
        >> $OUT_FILE

#Remove the temp file
    rm ${OUT_FILE}.temp

