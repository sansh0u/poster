import os
import pandas as pd
import re

def extract_gene_id(full_id):
    """提取基因 ID（例：3069at7088_1|... -> 3069at7088）"""
    if pd.isna(full_id):
        return None
    match = re.match(r'^([^_|]+at\d+)_', full_id)
    return match.group(1) if match else None

def process_luffia_comparisons():
    """仅处理 luffia__v__luffia1.tsv 与 luffia__v__Pachythelia.tsv，提取共享基因"""
    tsv_folder = os.path.join("orthofinder_input", "OrthoFinder", "Results_Apr28", "Orthologues", "Orthologues_luffia")

    # 固定两个文件路径
    ref_file = os.path.join(tsv_folder, "luffia__v__luffia1.tsv")
    other_file = os.path.join(tsv_folder, "luffia__v__Pachythelia.tsv")

    # 加载参考文件
    if not os.path.exists(ref_file) or not os.path.exists(other_file):
        print("缺少必要的 TSV 文件")
        return []

    df_ref = pd.read_csv(ref_file, sep='\t')
    df_other = pd.read_csv(other_file, sep='\t')

    df_ref['GeneID'] = df_ref['luffia'].apply(extract_gene_id)
    df_other['GeneID'] = df_other['luffia'].apply(extract_gene_id)

    ref_subset = df_ref[['Orthogroup', 'GeneID']].dropna()
    other_subset = df_other[['Orthogroup', 'GeneID']].dropna()

    merged = pd.merge(ref_subset, other_subset, on=['Orthogroup', 'GeneID']).drop_duplicates()

    return [{
        'species_file': "luffia__v__Pachythelia.tsv",
        'shared': merged
    }]

def extract_and_write_aln_files(shared_gene_data, busco_base_path="busco_outnew", output_folder="aln"):
    """从 BUSCO 结果中提取共享基因的 .faa 和 .fna 文件，重命名 header 并输出合并文件"""
    os.makedirs(output_folder, exist_ok=True)

    for entry in shared_gene_data:
        shared_df = entry['shared']

        for gene_id in shared_df['GeneID']:
            combined_faa = []
            combined_fna = []

            for species in ['luffia', 'luffia1', 'Pachythelia']:
                seq_folder = os.path.join(busco_base_path, species, "run_lepidoptera_odb10",
                                          "busco_sequences", "single_copy_busco_sequences")
                faa_path = os.path.join(seq_folder, f"{gene_id}.faa")
                fna_path = os.path.join(seq_folder, f"{gene_id}.fna")

                if os.path.exists(faa_path):
                    with open(faa_path) as f:
                        lines = f.read().splitlines()
                        if lines and lines[0].startswith(">"):
                            new_header = f">{species}_{gene_id}"
                            combined_faa.append("\n".join([new_header] + lines[1:]))

                if os.path.exists(fna_path):
                    with open(fna_path) as f:
                        lines = f.read().splitlines()
                        if lines and lines[0].startswith(">"):
                            new_header = f">{species}_{gene_id}"
                            combined_fna.append("\n".join([new_header] + lines[1:]))

            if len(combined_faa) == 3 and len(combined_fna) == 3:
                faa_out_path = os.path.join(output_folder, f"{gene_id}.faa")
                fna_out_path = os.path.join(output_folder, f"{gene_id}.fna")

                with open(faa_out_path, "w") as f:
                    f.write("\n".join(combined_faa) + "\n")

                with open(fna_out_path, "w") as f:
                    f.write("\n".join(combined_fna) + "\n")

if __name__ == "__main__":
    shared_gene_data = process_luffia_comparisons()
    if shared_gene_data:
        extract_and_write_aln_files(shared_gene_data)
