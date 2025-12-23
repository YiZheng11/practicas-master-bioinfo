#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH --output=logs/htseq_%j.log
#SBATCH --error=logs/htseq_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Exit if any command fails
set -e

# Añadir el bin de HTSeq al PATH?
export PATH=$HOME/.local/bin:$PATH

# Ver version de HTSeq
htseq-count --version

# Variables input
BAM_DIR=/home/proyectos/project_name/username/data/03_hisat2
GTF_FILE=/home/proyectos/project_name/username/data/GCF_000001635.27_GRCm39_MitoCarta.gtf
OUTPUT_DIR=/home/proyectos/project_name/username/04_htseq_filt_new
N_CORES=8

# First check if your env has GNU parallel ( UAM nodes should have it i guess)

if ! command -v parallel &> /dev/null; then
	echo "GNU parallel is not installed. Install it with: sudo apt install parallel"
	exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Lista de BAMs
BAM_FILES=($BAM_DIR/*_sorted.bam)

# Función para ejecutar HTSeq en un BAM
run_htseq() {
	BAM_FILE="$1"
	SAMPLE=$(basename "$BAM_FILE" _sorted.bam)
	OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE}.tsv"

	if [[ -s "$OUTPUT_FILE" ]]; then
		echo "Skipping $SAMPLE, it's already processed."
		continue
	fi

	echo "Processing $SAMPLE..."

	htseq-count \
		--format bam \
		--order pos \
		--mode intersection-strict \
		--stranded reverse \
		--minaqual 1 \
		--type gene \
		--idattr gene_id \
		"$BAM_FILE" "$GTF_FILE" > "$OUTPUT_FILE"
}

export -f run_htseq
export OUTPUT_DIR GTF_FILE

parallel -j $N_CORES run_htseq ::: "${BAM_FILES[@]}"
