set -xeu
set -o pipefail
indir=$1

# from docs
echo Obtaining Gold reference database for orienting detection
GOLD=https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip
if [ ! -e "${SYTYCS_CACHEDIR}/gold.fasta" ]
then

    if [ ! -e "${SYTYCS_CACHEDIR}/silva.gold.bacteria.zip" ]; then
        wget $GOLD -O ${SYTYCS_CACHEDIR}/silva.gold.bacteria.zip
    fi

    echo Decompressing and reformatting...
    unzip -p ${SYTYCS_CACHEDIR}/silva.gold.bacteria.zip silva.gold.align | \
        sed -e "s/[.-]//g" > ${SYTYCS_CACHEDIR}/gold.fasta

fi



cat $indir/pcr/*.fasta | sed "s| |--|g" > $indir/cluster_input.fasta
vsearch -orient $indir/cluster_input.fasta -db ${SYTYCS_CACHEDIR}/gold.fasta -fastaout $indir/cluster_input.oriented.fasta
vsearch --cluster_fast $indir/cluster_input.oriented.fasta --threads 16  --id 1 --strand both  -uc $indir/clusters.uc
mafft --auto $indir/cluster_input.fasta > $indir/alignment/amplicons.mafft.aln
# TODO: switch for infernal?

#FastTree -gtr -nt $indir/alignment/amplicons.mafft.aln > $indir/alignment/amplicons.mafft.aln.nwk
