
## A.*Pipeline de clasificación de archivos.gz por códigos de dengue en R*
### Establecer directorio de trabajo
setwd("/home/sec-dengue")

### Mostrar contenidos iniciales
cat("Archivos en la carpeta 'prueba':\n")
print(list.files("prueba"))

### Dentro de la carpeta pruebas, codigos.txts dengue
contenido1 <- readLines("prueba/dengue1.txt")
contenido2 <- readLines("prueba/dengue2.txt")
contenido3 <- readLines("prueba/dengue3.txt")
length(contenido1)
length(contenido2)
length(contenido3)

### Mostrar códigos cargados
cat("===== dengue1.txt =====\n")
print(contenido1)
cat("===== dengue2.txt =====\n")
print(contenido2)
cat("===== dengue3.txt =====\n")
print(contenido3)

### Crear carpetas 
dir.create("prueba/dengue1", showWarnings = FALSE)
dir.create("prueba/dengue2", showWarnings = FALSE)
dir.create("prueba/dengue3", showWarnings = FALSE)

### Buscar archivos .gz en carpeta prueba/
archivos <- list.files("prueba", pattern = "\\.gz$", full.names = TRUE)

### Contadores
cont1 <- 0
cont2 <- 0
cont3 <- 0

### copiar y separar cada archivo
for (archivo in archivos) {
  nombre <- basename(archivo)
  codigo <- substr(nombre, 1, 10)  # primeros 10 caracteres del nombre

  if (codigo %in% contenido1) {
    file.copy(archivo, file.path("prueba/dengue1", nombre), overwrite = TRUE)
    cont1 <- cont1 + 1
  } else if (codigo %in% contenido2) {
    file.copy(archivo, file.path("prueba/dengue2", nombre), overwrite = TRUE)
    cont2 <- cont2 + 1
  } else if (codigo %in% contenido3) {
    file.copy(archivo, file.path("prueba/dengue3", nombre), overwrite = TRUE)
    cont3 <- cont3 + 1
  }
}

cat("Archivos copiados a 'dengue1':", cont1, "\n")
cat("Archivos copiados a 'dengue2':", cont2, "\n")
cat("Archivos copiados a 'dengue3':", cont3, "\n")

## B.*Pipeline de Análisis de Lecturas NGS para Virus DENV2: Control de Calidad, Mapeo, Ensamblaje y Filtrado por Cobertura*

#!/bin/bash
# --- 1. Análisis de calidad con FastQC y MultiQC ---

mkdir -p fastqc_reporte

for file in *.fastq.gz; do
    echo "Procesando $file con FastQC..."
    xvfb-run fastqc "$file" -o fastqc_reporte/
done

cd fastqc_reporte
multiqc . -o .
xdg-open multiqc_report.html || echo "No se pudo abrir el reporte MultiQC automáticamente"
cd ..

# --- 2. Indexar el genoma de referencia ---

bwa index denv2_ref.fasta

# --- 3. Ensamblaje y procesamiento ---

for r1 in *L001_R1_001.fastq.gz; do
    prefix=$(basename "$r1" _L001_R1_001.fastq.gz)
    r2="${prefix}_L001_R2_001.fastq.gz"

    echo "Procesando muestra: $prefix"

    bwa mem -t 15 denv2_ref.fasta "$r1" "$r2" > "${prefix}_uno.sam"
    samtools view -@ 10 -bS -T denv2_ref.fasta "${prefix}_uno.sam" > "${prefix}_unoa.bam"
    samtools sort -@ 10 -n "${prefix}_unoa.bam" -o "${prefix}_dosa.bam"
    samtools fixmate -@ 10 -m "${prefix}_dosa.bam" "${prefix}_tresa.bam"
    samtools sort -@ 10 "${prefix}_tresa.bam" -o "${prefix}_cuatroa.bam"
    samtools markdup -@ 10 "${prefix}_cuatroa.bam" "${prefix}.bam"
    samtools index -@ 10 "${prefix}.bam"

    rm "${prefix}_uno.sam" "${prefix}_unoa.bam" "${prefix}_dosa.bam" "${prefix}_tresa.bam" "${prefix}_cuatroa.bam"
done

# --- 4. Obtener consenso y variantes con código corto (10 caracteres) ---

for bam in *.bam; do
    prefix=$(basename "$bam" .bam)
    codigo=${prefix:0:10}

    echo "Generando consenso y variantes para: $codigo"

    samtools mpileup -aa -A -d 0 -Q 0 "$bam" \
    | ivar consensus -p "$codigo" -q 25 -t 0.6 -m 10

    samtools mpileup -aa -A -d 600000 -B -Q 0 "$bam" \
    | ivar variants -p "${codigo}.tsv" -q 25 -t 0.6 -r denv2_ref.fasta
done

# --- 5. Corregir encabezados y renombrar archivos .fa con código corto ---

for fa in *.fa; do
    base=$(basename "$fa" .fa)
    codigo=${base:0:10}
    echo "Corrigiendo encabezado en: $fa con código: $codigo"
    sed -i "1s/.*/>${codigo}/" "$fa"

    if [[ "$fa" != "${codigo}.fa" ]]; then
        if [[ -f "${codigo}.fa" ]]; then
            echo "Advertencia: ${codigo}.fa ya existe. No se renombró $fa"
        else
            mv "$fa" "${codigo}.fa"
        fi
    fi
done

# Combinar todos los .fa en un multifasta
cat *.fa > secuencias.fasta

# --- 6. Estimar profundidad y cobertura con código corto ---

echo -e "Codigo\tProfundidad(x)\tCobertura(%)" > profundidad_ns.txt

for bam in *.bam; do
    prefix=$(basename "$bam" .bam)
    codigo=${prefix:0:10}

    DEP=$(samtools depth "$bam" | awk '{sum+=$3} END {printf "%.2f", sum/10723}')

    if [[ -f "${codigo}.fa" ]]; then
        COV=$(seqtk comp "${codigo}.fa" | awk '{x=$3+$4+$5+$6; y=10723; printf "%.2f", (1-(y-x)/y)*100}')
    else
        COV="0.00"
    fi

    echo -e "${codigo}\t${DEP}x\t${COV}%" >> profundidad_ns.txt
done

# --- 7. Filtrar secuencias con cobertura > 90% y extraer del multifasta ---

# 1. Extraer códigos con cobertura mayor a 90%
awk -F'\t' 'NR>1 && $3+0 > 90 { print $1 }' profundidad_ns.txt > secuencias_mayores90.txt

# 2. Contar cuántas cumplen el criterio
TOTAL=$(wc -l < secuencias_mayores90.txt)
echo " Se encontraron $TOTAL secuencias con cobertura mayor al 90%."

# 3. Extraer esas secuencias del archivo multifasta
seqtk subseq secuencias.fasta secuencias_mayores90.txt > secuencias_mayores90.fasta
grep -c '^>' secuencias_mayores90.fasta
done































