#!/bin/bash
#SBATCH --cpus-per-task=4

# Añadir Trim Galore al PATH (ajusta la ruta a donde lo instalaste)
export PATH=$HOME/.local/bin:$PATH

module load cutadapt/4.1

# Comprobar Trim Galore
trim_galore --version

# Directorios
RAW_DIR=/home/proyectos/project_name/username/data/00_reads
OUT_DIR=/home/proyectos/project_name/username/01_trimmed

# Loop sobre todos los pares de fastq
for R1 in "$RAW_DIR"/*9a_1.fq.gz; do
    BASE=$(basename "$R1" "_1.fq.gz")
    R2="$RAW_DIR/${BASE}_2.fq.gz"

    # Archivos de salida
    OUT_R1="$OUT_DIR/${BASE}_1_val_1.fq.gz"
    OUT_R2="$OUT_DIR/${BASE}_2_val_2.fq.gz"

    # Saltar si ya existen
    if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
        echo "Skipping $BASE, it already exists."
        continue
    fi

    echo "Processing $BASE..."

    trim_galore --paired \
		--cores 4 \
		-o "$OUT_DIR" \
		$R1 $R2
done
