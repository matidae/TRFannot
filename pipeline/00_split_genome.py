import subprocess
from pathlib import Path
import yaml

# python pipeline/00_split_genome.py config.yaml
cfg      = yaml.safe_load(open("config.yaml"))
workdir  = Path(cfg["workdir"])
fasta    = Path(cfg["genome"]["fasta"])
fai      = Path(f"{fasta}.fai")
chunk_size = cfg["chunking"]["chunk_size"]
min_length = cfg["chunking"]["min_length"]
chunks_dir = workdir / cfg["chunking"]["out_dir"]

chunks_dir.mkdir(parents=True, exist_ok=True)

# Parse fai and yield (chrom, start, end) for sequences above min_length
def make_chunks(fai, chunk_size, min_length):
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    for line in open(fai):
        chrom, length = line.split("\t")[:2]
        length = int(length)
        if length < min_length:
            continue
        start = 0
        while start < length:
            end = min(start + chunk_size, length)
            yield chrom, start, end
            start = end

# samtools faidx data/genome/LisVul1.1.fas CHR1:1-100000000 > data/chunks/CHR1.1-100000000.fa
def extract_chunk(fasta, chrom, start, end, out_path):
    region = f"{chrom}:{start + 1}-{end}"
    with open(out_path, "w") as f:
        result = subprocess.run(["samtools", "faidx", str(fasta), region], stdout=f)
    if result.returncode != 0:
        print(f"ERROR: samtools failed on {region}")


chunks = list(make_chunks(fai, chunk_size))
print(f"Extracting {len(chunks)} chunks")


for chrom, start, end in chunks:
    out_path = chunks_dir / f"{chrom}.{start + 1}-{end}.fa"
    extract_chunk(fasta, chrom, start, end, out_path)
    print(f"  {out_path.name}")
print("Splitting genome done")