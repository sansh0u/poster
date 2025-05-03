#!/bin/bash

cd aln

for aln in *.aln; do
    base="${aln%.aln}"
    fna="${base}.fna"
    out="${base}.fasta"

    if [[ -f "$fna" ]]; then
        echo "正在处理 $base..."
        perl pal2nal.pl "$aln" "$fna" -output fasta > "$out"
    else
        echo "❌ 缺少 $fna，跳过 $aln"
    fi
done

