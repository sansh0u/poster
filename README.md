# **1.Data**

 All genome files were obtained from NCBI.

# **2.Method**

 1. Annotating genomes using BUSCO (version 5.8.2). Using the lepidoptera_odb10 database and using the metaeuk pipeline.

```bash
 busco --in genome.fna --out busco_outnew --lineage_dataset lepidoptera_odb10 --mode genome --metaeuk --cpu 16
```

2. In the busco result, extract the single copy gene of each species, and place them in a folder. Using OrthoFinder (version 3.0.1b1) to find single copy orthologous genes.

```bash
orthofinder -f your_data_folder
```
3. Using MAFFT (version 7.526) to align protein sequences, and then using pal2nal to map the alignment back to nucleotide sequences. The Shell script file "run_mafft.sh" and "run_pal2nal.sh" were used for batch operations.

4. Using self-made py files for McDonald–Kreitman Test.

```bash
python test1.py --dir your_dircrtion --out mkt.tsv
```
5. Using Emapper(version 2.1.12) and eggNOG database (version 5.0.2) to perform GO functional annotation on all orthologous genes of luffia.

```bash
 emapper.py -i combined.faa -o luffia_pos --itype proteins --data_dir ~/eggnog_data --cpu 16
 ```
   
6. Using RStudio to filter the mkt results. Delete all genes without GO function，Positively selected genes were filtered based on the criteria: α > 0, NI < 1, and Fisher test p-value < 0.05. Negative selected genes were filtered based on the criteria: α < 0, NI > 1, and Fisher test p-value < 0.05.

7. The positive selection gene was used as the target gene, and all orthogonal genes were used as background genes. RStudio was used to perform GO term analysis. The negative selection gene was analyzed in the same way.
