#Downloads and converts sigProfiler (COSMIC v3.1) to TSV.
#Other sigProfiler signatures need to be manually downloaded from Synapse
#and manipulated to the required format

#File name variables
    IN_SHEET="https://cancer.sanger.ac.uk/sigs-assets-20/COSMIC_Mutational_Signatures_v3.1.xlsx"

    OUT_DIR="data/raw/reference_signatures/"
    OUT_PREFIX="sigProfiler_"
    OUT_SUFFIX="_signatures.tsv"

#Make sure the output directory exists
    mkdir -p $OUT_DIR

#Download
    wget -nv $IN_SHEET -O ${OUT_DIR}.temp.cosmic.xlsx

#COSMIC v3
    #Read sheets to CSVs
    ssconvert ${OUT_DIR}.temp.cosmic.xlsx -S ${OUT_DIR}.temp.cosmic.csv

    #SBS
    #Sort data by mutation, then context and
    #replace column name Subtype -> SubType and replace commas with tabs
    cat <(head -1 ${OUT_DIR}.temp.cosmic.csv.0) \
        <(tail -n +2 ${OUT_DIR}.temp.cosmic.csv.0 | sort) |
        perl -pe 's/(^|,)Subtype(,|$)/$1SubType$2/; s/,/\t/g' \
        > ${OUT_DIR}${OUT_PREFIX}SBS${OUT_SUFFIX}

    #DBS
    #Replace commas with tabs
    perl -pe 's/,/\t/g' ${OUT_DIR}.temp.cosmic.csv.1 \
        > ${OUT_DIR}${OUT_PREFIX}DBS${OUT_SUFFIX}

    #ID
    #Remove unused types
    grep -vP "^(\d:Ins:M:\d|complex|non_matching)," ${OUT_DIR}.temp.cosmic.csv.2 |
        #Change the type format
        perl -ne '
            if ($. > 1) {
                s/^([^,]+),//;
                @F = split /:/, $1;

                $type = uc($F[1]) . "_";

                if ($F[2] =~ m/^[CT]$/) {
                    $type .= "$F[2]";
                } elsif ($F[2] eq "R") {
                    $type .= "repeats";
                } elsif ($F[2] eq "M") {
                    $type .= "MH";
                } else {
                    die "Unexpected mutation type " . join(":", @F) . "!\n"
                }

                $type .= "_$F[0]_$F[3]";

                $type =~ s/5/5\+/g;

                $_ = "$type,$_";
            }
            print;
        ' |
        #Replace commas with tabs
        perl -pe 's/,/\t/g' \
            > ${OUT_DIR}${OUT_PREFIX}ID${OUT_SUFFIX}

#Remove the temp files
    rm ${OUT_DIR}.temp.cosmic.xlsx ${OUT_DIR}.temp.cosmic.csv.*

