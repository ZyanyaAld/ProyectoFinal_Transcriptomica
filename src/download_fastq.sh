#!/bin/bash
#SBATCH --job-name=downloadRNA
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00

cd /export/space3/users/zyanyava/2026-2/Transcriptomica/RNA_seq_Project/data/raw

URLS=(
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/068/SRR25629468/SRR25629468_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/068/SRR25629468/SRR25629468_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/073/SRR25629473/SRR25629473_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/073/SRR25629473/SRR25629473_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/067/SRR25629467/SRR25629467_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/067/SRR25629467/SRR25629467_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/070/SRR25629470/SRR25629470_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/070/SRR25629470/SRR25629470_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/071/SRR25629471/SRR25629471_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/071/SRR25629471/SRR25629471_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/065/SRR25629465/SRR25629465_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/065/SRR25629465/SRR25629465_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/066/SRR25629466/SRR25629466_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/066/SRR25629466/SRR25629466_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/069/SRR25629469/SRR25629469_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/069/SRR25629469/SRR25629469_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/072/SRR25629472/SRR25629472_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR256/072/SRR25629472/SRR25629472_2.fastq.gz
)

download() {

    local url=$1
    local file=$(basename "$url")

    if [[ -f "$file" ]]; then
        echo "[SKIP] $file already exists"
        return
    fi

    echo "[DOWNLOADING] $file"

    wget -c "$url"

    if [[ $? -eq 0 ]]; then
        echo "[SUCCESS] $file"
    else
        echo "[FAILED] $file"
    fi
}

export -f download

printf '%s\n' "${URLS[@]}" | xargs -P 2 -I {} bash -c 'download "$@"' _ {}

echo "Download complete!"