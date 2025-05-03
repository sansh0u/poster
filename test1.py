import os
import argparse
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
from itertools import combinations
from typing import List, Dict
from scipy.stats import fisher_exact
# mkt_focal.py：以 focal 物种为主物种，分别与内群 & 外群计算 Pn/Ps 与 Dn/Ds

# 翻译密码子为氨基酸
def translate_codon(codon: str) -> str:
    try:
        return str(Seq(codon).translate(table=1))
    except:
        return "X"

# 判断是否为同义突变
def is_synonymous(codon1: str, codon2: str) -> bool:
    return translate_codon(codon1) == translate_codon(codon2) and 'N' not in codon1+codon2 and '-' not in codon1+codon2

# 从序列名中解析各物种索引
def parse_species(names: List[str], main: str, ingroup: str, outgroup: str):
    main_idx = []
    ingroup_idx = []
    outgroup_idx = []
    for i, name in enumerate(names):
        prefix = name.split('_')[0]
        if prefix == main:
            main_idx.append(i)
        elif prefix == ingroup:
            ingroup_idx.append(i)
        elif prefix == outgroup:
            outgroup_idx.append(i)
    return main_idx, ingroup_idx, outgroup_idx

# 执行 focal MKT 计算
def focal_mkt(aln, main_idx, ingroup_idx, outgroup_idx) -> Dict[str, float]:
    Pn = Ps = Dn = Ds = skipped = 0
    aln_len = aln.get_alignment_length()
    num_codons = aln_len // 3

    for i in range(num_codons):
        codon_pos = slice(i * 3, (i + 1) * 3)

        # 提取单个密码子（注意每个组只有一个样本）
        main_codon = str(aln[main_idx[0], codon_pos].seq).upper()
        in_codon = str(aln[ingroup_idx[0], codon_pos].seq).upper()
        out_codon = str(aln[outgroup_idx[0], codon_pos].seq).upper()

        all_codons = main_codon + in_codon + out_codon
        if 'N' in all_codons or '-' in all_codons:
            skipped += 1
            continue

        # 群内多态性：main vs ingroup
        if main_codon != in_codon:
            if is_synonymous(main_codon, in_codon):
                Ps += 1
            else:
                Pn += 1

        # 固定差异：main vs outgroup（要求两边一致内部无变异且不同）
        if main_codon != out_codon:
            if is_synonymous(main_codon, out_codon):
                Ds += 1
            else:
                Dn += 1

    result = {'Pn': Pn, 'Ps': Ps, 'Dn': Dn, 'Ds': Ds, 'skipped_codons': skipped}
    if Ps > 0 and Ds > 0:
        result['NI'] = (Pn / Ps) / (Dn / Ds)
        result['alpha'] = 1 - result['NI']
    else:
        result['NI'] = None
        result['alpha'] = None
    
    if all(x >= 1 for x in [Pn, Ps, Dn, Ds]):
        try:
            table = [[Dn, Ds], [Pn, Ps]]
            _, pval = fisher_exact(table, alternative='two-sided')
            result['pvalue'] = pval
        except:
            result['pvalue'] = None
    else:
        result['pvalue'] = None

    return result

def main():
    parser = argparse.ArgumentParser(description="Focal-species MKT analysis.")
    parser.add_argument("--fasta", required=True, help="Codon-aligned FASTA file")
    parser.add_argument("--main", required=True, help="Prefix of focal species (e.g. luffia)")
    parser.add_argument("--ingroup", required=True, help="Prefix of sister species for polymorphism (e.g. luffia1)")
    parser.add_argument("--outgroup", required=True, help="Prefix of outgroup species (e.g. Pachythelia)")
    args = parser.parse_args()

    aln = AlignIO.read(args.fasta, "fasta")
    names = [record.id for record in aln]
    main_idx, ingroup_idx, outgroup_idx = parse_species(names, args.main, args.ingroup, args.outgroup)

    print("Focal species indices:", main_idx)
    print("Sister group indices:", ingroup_idx)
    print("Outgroup indices:", outgroup_idx)

    if not main_idx or not ingroup_idx or not outgroup_idx:
        print("Error: one or more groups not found.")
        return

    stats = focal_mkt(aln, main_idx, ingroup_idx, outgroup_idx)
    print("MKT (focal) results:")
    for k, v in stats.items():
        print(f"{k}: {v}")


def batch_mkt_focal(directory: str, main: str, ingroup: str, outgroup: str, suffix: str = ".fasta") -> pd.DataFrame:
    results = []
    for fname in os.listdir(directory):
        if not fname.endswith(suffix):
            continue
        path = os.path.join(directory, fname)
        try:
            aln = AlignIO.read(path, "fasta")
            names = [record.id for record in aln]
            main_idx, ing_idx, out_idx = parse_species(names, main, ingroup, outgroup)

            if not main_idx or not ing_idx or not out_idx:
                print(f"Skipping {fname}: missing group.")
                continue

            stats = focal_mkt(aln, main_idx, ing_idx, out_idx)
            stats["file"] = fname
            results.append(stats)
        except Exception as e:
            print(f"Error processing {fname}: {e}")

    return pd.DataFrame(results)


# 命令行接口
def main():
    parser = argparse.ArgumentParser(description="Batch focal MKT analysis.")
    parser.add_argument("--dir", required=True, help="Directory with fasta files")
    parser.add_argument("--main", required=True, help="Prefix of focal species")
    parser.add_argument("--ingroup", required=True, help="Prefix of sister species (ingroup)")
    parser.add_argument("--outgroup", required=True, help="Prefix of outgroup species")
    parser.add_argument("--out", required=True, help="Output CSV file path")
    args = parser.parse_args()

    df = batch_mkt_focal(args.dir, args.main, args.ingroup, args.outgroup)
    df.to_csv(args.out, index=False)
    print(f"Results written to {args.out}")
if __name__ == "__main__":
    main()
