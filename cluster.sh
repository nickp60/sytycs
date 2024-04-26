set -xeu
set -o pipefail
indir=$1
cat $indir/pcr/*.fasta > $indir/pcr/input.fasta
vsearch --cluster_fast $indir/pcr/input.fasta --threads 16  --id 1 --strand both  -uc $indir/pcr/input.fasta.uc
