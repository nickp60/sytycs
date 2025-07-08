#!/usr/bin/env python
import os
import sys
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def parse_primersearch(f):
    """ This splits apart records from primersearch; currently being done with awk

    """
    pass


def main():
    fasta, psout, maxlen = sys.argv[1], sys.argv[2], sys.argv[3]
    with gzip.open(fasta, "rt") as handle:
        records_pre = list(SeqIO.parse(handle, "fasta"))
        recids = [x.id  for x in records_pre]
        records = {k:v for k,v in zip(recids, records_pre)}
#        recdict =
    with open(psout, "r") as inf:
        for line in inf:
            if line == "\n":
                break
            sys.stderr.write(",".join(line.split("\t")) )
            amplimer, acc, description, start, revstop, amplimer_len = line.strip().split("\t")
            if int(amplimer_len) > maxlen:
                sys.stderr.write(f"warning: amplimer {amplimer} exceeds {maxlen} bp; depending on your primers this might be a spurious hit.`\n")
                continue
            acc = acc.replace("Sequence: ", "").split()[0]
            rec = records[acc]
            true_stop = len(rec.seq) - int(revstop.replace("[", "").replace("]", ""))
            name = f"{rec.id}:{int(start)}-{true_stop}"
            amplicon = rec.seq[int(start):true_stop]
            if amplicon.startswith("CGTCA"):
                amplicon = amplicon.reverse_complement()
                name = name + "-RC"

            sys.stdout.write(SeqRecord(id=name, description=description, seq=amplicon).format("fasta") )

if __name__ == "__main__":
    main()
