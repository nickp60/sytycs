set -xeu
set -o pipefail
indir=$1

# from docs
echo Obtaining Gold reference database for orienting detection
GOLD=https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip
if [ ! -f "${SYTYCS_CACHEDIR}/gold.fasta" ]
then

    if [ ! -f "${SYTYCS_CACHEDIR}/silva.gold.bacteria.zip" ]; then
        wget $GOLD -O ${SYTYCS_CACHEDIR}/silva.gold.bacteria.zip
    fi

    echo Decompressing and reformatting...
    unzip -p ${SYTYCS_CACHEDIR}/silva.gold.bacteria.zip silva.gold.align | \
        sed -e "s/[.-]//g" > ${SYTYCS_CACHEDIR}/gold.fasta

fi



cat $indir/pcr/*.fasta | sed "s| |--|g" > $indir/cluster_input.fasta
#vsearch --makeudb_usearch ${SYTYCS_CACHEDIR}/gold.fasta --output ${SYTYCS_CACHEDIR}/gold.udb
#vsearch -orient $indir/cluster_input.fasta -db ${SYTYCS_CACHEDIR}/gold.udb -fastaout $indir/cluster_input.oriented.fasta

vsearch -orient $indir/cluster_input.fasta -db ${SYTYCS_CACHEDIR}/gold.fasta -fastaout $indir/cluster_input.oriented.fasta

vsearch --cluster_fast $indir/cluster_input.oriented.fasta --threads 16  --id 1 --strand both  -uc $indir/clusters.uc

# add in E. coli to root this
echo -n ">E. coli\nACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGG" >>  $indir/cluster_input.oriented.fasta
mafft --auto $indir/cluster_input.oriented.fasta > $indir/alignment/amplicons.mafft.aln
# TODO: switch for infernal?
# the sed prevents errors where fasttree only takes the sequence id prior to the :, leading to duplicates
cat $indir/alignment/amplicons.mafft.aln | sed "s|:|-|g" | FastTree -gtr -nt  > $indir/alignment/amplicons.mafft.aln.nwk
