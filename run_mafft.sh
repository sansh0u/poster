#!/bin/bash
cd aln

for file in *.faa; do
    out="${file%.faa}.aln"
    echo "Aligning $file -> $out"
    mafft --maxiterate 1000 --globalpair "$file" > "$out"
done
