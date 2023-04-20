#Prepares data for annotating transcription strand on mutations.

#Makes GENCODE v29 annotation bed file with only protein coding gene and exon
#rows and with duplicates and info columns are removed. Then creates BED files
#of transcripts and of exons with strand information. Finally creates a table
#of chromosome sizes.

#File name variables
    IN_FILE="data/raw/genome_annotations/GRCh38.GENCODE.v29.gtf"
    REF_FASTA="data/raw/reference_genomes/GRCh38.d1.vd1.fa"

    GENES_DIR="data/interim/genome_annotations/"
    GENES_FILE="${GENES_DIR}GRCh38.GENCODE.v29.known_genes.bed"
    GENES_TRANSCRIPTS="data/interim/genome_annotations/GRCh38.GENCODE.v29.known_genes_whole.bed"
    GENES_EXONS="data/interim/genome_annotations/GRCh38.GENCODE.v29.known_genes_exons.bed"
    CHR_DIR="data/interim/reference_genomes/"
    CHR_SIZES="${CHR_DIR}GRCh38.chr_sizes.txt"
    CHR_FILE="${CHR_DIR}GRCh38.chr_order.txt"

#NEED TO DECIDE HOW EXACTLY TO FILTER / CHOOSE REGIONS

#ALEXANDROV ET AL. 2013:
#"Strand bias catalogs were derived for each sample using only substitutions
#identified in the transcribed regions of well-annotated protein coding genes.
#Genomic regions of bidirectional transcription were excluded from the strand
#bias analysis."

#ALEXANDROV ET AL. 2020:
#"-- mutations within transcribed genome regions were selected and classified
#according to whether the pyrimidine of the mutated base pair fell on the
#transcribed or untranscribed strand --."
#Based on the GENCODE BED file in preprocessing code folders, "transcripts" of
#GENCODE v19/hs37 or v20/GRCh38 was used.

#Make sure output directories exist
    mkdir -p $GENES_DIR $CHR_DIR

#Get desired entries with relevant information only and output as BED
    grep -v "^##"  $IN_FILE |
        #Take any transcripts or exons (all are nested in transcripts)
        grep -P "\t(transcript|exon)\t" |
        #Keep only genes/exons (rows) of well annotated protein genes (levels 1,2)
        #grep -P "\t(gene|exon)\t.*(gene_type \"protein_coding\"; gene_status \"KNOWN\";).*level [1-2];" |
        #Remove unnecessary columns, incl. score and information columns
        perl -ane '
            ($gene_name) = $_ =~ /gene_name\ \"([^\"]*)\"\;/;
            @output = (@F[(0,3,4)], $gene_name, @F[(6,2)]);
            print(join("\t", @output) . "\n");
        ' |
        #Sort entries
        sed 's/^chrM/chrZ/' |
        sort -k1,1V -k2,2n -k3,3n |
        sed 's/^chrZ/chrM/' |
        #Take exon only once (remove duplicates from different transcripts)
        uniq \
        > $GENES_FILE &

#Create chromosome sizes file
    grep "^>" $REF_FASTA |
        perl -pe 's/^>//; s/\s+.*LN\:(\d+).*$/\t$1/' \
        > $CHR_SIZES

wait

#Create whole transcripts / all exons split bed files with intervals merged by
#strand
    #Whole transcripts
    grep "transcript" $GENES_FILE |
        #Make sure strand is sixth column
        awk '{temp=$5; $5=$6; $6=temp; print}' OFS='\t' |
        #Merge contiguous/overlapping intervals on the same strand
        bedtools merge -i stdin -s -c 4,5,6 -o distinct |
        #Sort entries
        sed 's/^chrM/chrZ/' |
        sort -k1,1V -k2,2n -k3,3n |
        sed 's/^chrZ/chrM/' \
        > $GENES_TRANSCRIPTS
    #All exons
    grep "exon" $GENES_FILE |
        #Remove transcripts with 'exon' in the gene name
        grep -v "transcript" |
        #Make sure strand is sixth column
        awk '{temp=$5; $5=$6; $6=temp; print}' OFS='\t' |
        #Merge contiguous/overlapping intervals on the same strand
        bedtools merge -i stdin -s -c 4,5,6 -o distinct |
        #Sort entries
        sed 's/^chrM/chrZ/' |
        sort -k1,1V -k2,2n -k3,3n |
        sed 's/^chrZ/chrM/' \
        > $GENES_EXONS

#Create chromosome order file with chromosome sizes
    #Empty chromosome list file
    rm -f $CHR_FILE

    #Create chromosome list with sizes
    while read CHR; do
        grep -P "$CHR\t" $CHR_SIZES \
            >> $CHR_FILE
    done < <( cut -f1 $GENES_FILE | uniq )

