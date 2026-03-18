#!/usr/bin/env python3

import glob
import os
import re
from collections import defaultdict

# Offsets for split chromosomes (from samtools faidx)
CHR_OFFSETS = {
    "CHR_1_1": 0,
    "CHR_1_2": 2144885605,

    "CHR_2_1": 0,
    "CHR_2_2": 2136077662,

    "CHR_3_1": 0,
    "CHR_3_2": 2137795666,

    "CHR_4_1": 0,
    "CHR_4_2": 2145962954,
}

def normalize_chrom(chrom):
    """
    Return (biological_chrom, offset)
    """
    if chrom in CHR_OFFSETS:
        base = chrom.rsplit("_", 1)[0]   # CHR_1_2 -> CHR_1
        return base, CHR_OFFSETS[chrom]
    else:
        return chrom, 0


def main():
    records_by_chr = defaultdict(list)

    for fname in glob.glob("trf_out/*.coords.tsv"):
        with open(fname) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue

                fields = line.split("\t")
                chrom = fields[0]
                start = int(fields[1])
                end   = int(fields[2])

                bio_chrom, offset = normalize_chrom(chrom)

                # Recalculate coordinates if needed
                start += offset
                end   += offset

                fields[0] = bio_chrom
                fields[1] = str(start)
                fields[2] = str(end)

                records_by_chr[bio_chrom].append(fields)

    # Write one file per chromosome
    os.makedirs("trf_merged", exist_ok=True)

    for chrom, records in records_by_chr.items():
        outname = f"trf_mergend/{chrom}.trf.tsv"
        with open(outname, "w") as out:
            for fields in records:
                out.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
