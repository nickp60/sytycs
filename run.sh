set -eux
set -o pipefail
if [ ! -f prokaryotes.txt ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
fi

## subset to "Complete Genome" for now.  Could include Chromosome? defined as https://www.ncbi.nlm.nih.gov/datasets/docs/v2/glossary/
## Chromosome—There is a sequence for one or more chromosomes. This may be a completely sequenced chromosome without gaps or a chromosome containing scaffolds or contigs with gaps between them. There may also be unplaced or unlocalized scaffolds.
## Complete genome—All chromosomes are gapless and contain runs of nine or less ambiguous bases (Ns), there are no unplaced or unlocalized scaffolds, and all the expected chromosomes are present (i.e., the assembly is not noted as having partial genome representation). Plasmids and organelles may or may not be included in the assembly, but if they are present, the sequences are gapless.
query=$1
primerf=$2
primerr=$3
outdir=sytycs_${query}-$primerf-$primerr/
mkdir -p $outdir
mkdir -p $SYTYCS_CACHEDIR
mkdir -p $outdir/pcr/
mkdir -p  $outdir/alignment/


# subset to lines where the first column starts with the query
# and to where the ftp field is not "-"
cat prokaryotes.txt | awk -F "\t"  '{if ($16 == "Complete Genome")  { print }}' | awk -F "\t"  '{if ($21 != "-")  { print }}' | awk -v query=$query -F "\t" '$1 ~ "^" query' > $outdir/tmp.queryhits.tsv
echo "      ------   found $(wc -l tmp.queryhits.tsv) hits for $query ------- "
echo "      ------   downloading  genomes ------- "
cat $outdir/tmp.queryhits.tsv | cut -f21 | while read urlpre; do
    url=${urlpre}_genomic.fna.gz
    base=$(basename $url)
    target="$SYTYCS_CACHEDIR/$base"
    if [ ! -f "$target" ]
    then
	# try to download it, but in cases like GCF_002278015.2_ASM227801v2 there may be multiple versions organised as subdirectories
	assemblybase=$(basename ${urlpre})
        alturl=${urlpre}/${assemblybase}_genomic.fna.gz
      wget   $url -O $target || wget $alturl -O $target
    else
	echo "                  ---  already have $target      "
    fi
    echo $target >> $outdir/manifest.txt
done

## create the primersearch input file:
echo -e "query\t$primerf\t$primerr" > $outdir/primers.txt
if [ -f "${outdir}pcr/tmp.primersearch" ]
then
rm ${outdir}pcr/tmp.primersearch
fi

echo "      ------   running primersearch on  genomes ------- "
cat $outdir/manifest.txt | while read infile
do
    if [ -s "$infile" ]
    then
	outfile=${outdir}pcr/$(basename $infile).fasta
    zcat $infile | seqret -filter -osformat embl >> tmp.embl
    # awk '{$1=$1;print}' trims whitespace through awk voodoo
    # | sed "s|.* strand at ||g" trims off the match details except for the corrdinate
    # the sed calls trim up the awful format to something tabular
    primersearch -infile $outdir/primers.txt -seqall tmp.embl -outfile ${outfile}.ps -mismatchpercent 5
    awk '{$1=$1;print}' ${outfile}.ps | tail -n+3 | sed "s|.* strand at ||g" | tr "\n" "\t"  | \
        awk 'BEGIN{RS="bp"}{print}' | sed "s|^\t||g" | sed "s|Amplimer length: ||g" | sed "s|with . mismatches||g" > ${outfile}.hits

    #this filters for hits shorter than 400, trims off brackets and extra whitespace, extracts the forward and reverse coords, separates them with ..,
    # and turns new lines to commas to make a comma-separated list
    # the reverse coords are, inconveniently, relative to the reverse strand :(
    # I wish I was joking...
#    coord_string=$(cat ${outfile}.hits | awk -F "\t" '{if ($6 <= 400) {print}}'  | grep Sequence |  tr -d "]" | tr -d "["  | cut -f 4,5 | sed "s|\t|..|g" | tr -d " " | tr "\n" ",")
#    coord_string=$(cat ${outfile}.hits | awk -F "\t" '{if ($6 <= 400) {print}}'  | grep Sequence |  awk 'BEGIN{FS=OFS="\t"} {print $4, $6+$4}'  | sed "s|\t|..|g" | tr -d " " | tr "\n" ",")
#    extractseq tmp.embl -reg "$coord_string" stdout -separate  > $outfile

    ./primersearch_to_fasta.py $infile  ${outfile}.hits >  ${outfile}
    rm tmp.embl
    else
echo "$infile doesn't exist or is empty"
fi

done
