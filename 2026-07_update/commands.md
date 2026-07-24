Run the assemblies:
```bash
threads=32

declare -A genome_sizes
genome_sizes["Enterobacter_hormaechei"]=5384747
genome_sizes["Klebsiella_pneumoniae"]=5990196
genome_sizes["Listeria_innocua"]=2972545
genome_sizes["Providencia_rettgeri"]=4465806
genome_sizes["Shigella_flexneri"]=4828487

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    genome_size=$genome_sizes["$s"]
    for i in {1..6}; do
        /usr/bin/time -v -o assemblies/myloasm_v0.6.0_"$i".time autocycler helper myloasm --reads reads_subsampled/"$i".fastq --out_prefix assemblies/myloasm_v0.6.0_"$i" --threads "$threads" --genome_size "$genome_size" --min_depth_rel 0.1
        /usr/bin/time -v -o assemblies/metamdbg_v1.4_"$i".time autocycler helper metamdbg --reads reads_subsampled/"$i".fastq --out_prefix assemblies/metamdbg_v1.4_"$i" --threads "$threads" --genome_size "$genome_size" --min_depth_rel 0.1
        /usr/bin/time -v -o assemblies/ilesta_v1.2.1_"$i".time autocycler helper ilesta --reads reads_subsampled/"$i".fastq --out_prefix assemblies/ilesta_v1.2.1_"$i" --threads "$threads" --genome_size "$genome_size" --min_depth_rel 0.1
    done
done

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    genome_size=$genome_sizes["$s"]
    for i in {1..6}; do
        /usr/bin/time -v -o assemblies/hybracter_v0.14.0_"$i".time hybracter long-single -l reads_subsampled/"$i".fastq -s "$s"_"$i" --auto -o hybracter_v0.14.0_"$i" -t 32 --skip_qc
        cp hybracter_v0.14.0_"$i"/FINAL_OUTPUT/*complete/*_final.fasta assemblies/hybracter_v0.14.0_"$i".fasta
        rm -r hybracter_v0.14.0_"$i"
    done
done
```

Run Autocycler pipelines:
```bash
threads=32

declare -A genome_sizes
genome_sizes["Enterobacter_hormaechei"]=5384747
genome_sizes["Klebsiella_pneumoniae"]=5990196
genome_sizes["Listeria_innocua"]=2972545
genome_sizes["Providencia_rettgeri"]=4465806
genome_sizes["Shigella_flexneri"]=4828487

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for i in {1..6}; do
        cd ~/2025-04_Autocycler_paper/"$s"

        mkdir autocycler_"$i"_2026-07 && cd autocycler_"$i"_2026-07
        /usr/bin/time -v -o ../assemblies/autocycler_2026-07_"$i".time autocycler_full.sh ../reads_subsampled/"$i".fastq 8 4
        cp autocycler_out/consensus_assembly.fasta ../assemblies/autocycler_2026-07_"$i".fasta

        mkdir autocycler_fast_"$i"_2026-07 && cd autocycler_fast_"$i"_2026-07
        /usr/bin/time -v -o ../assemblies/autocycler_fast_2026-07_"$i".time autocycler_fast.sh ../reads_subsampled/"$i".fastq 8 4
        cp autocycler_out/consensus_assembly.fasta ../assemblies/autocycler_fast_2026-07_"$i".fasta
    done
done
```

Create a results file for these new assemblies:
```bash
cd ~/2025-04_Autocycler_paper
./assess_assembly.py --header > results_asm5.tsv
./assess_assembly.py --header > results_asm10.tsv
./assess_assembly.py --header > results_asm20.tsv
./assess_assembly.py --header > results_map-ont.tsv

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/myloasm_v0.6.0*.fasta "$s"/assemblies/metamdbg_v1.4*.fasta "$s"/assemblies/ilesta_v1.2.1*.fasta "$s"/assemblies/hybracter_v0.14.0*.fasta "$s"/assemblies/autocycler_fast_2026-07*.fasta "$s"/assemblies/autocycler_2026-07*.fasta; do
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset asm5 >> results_asm5.tsv
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset asm10 >> results_asm10.tsv
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset asm20 >> results_asm20.tsv
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset map-ont >> results_map-ont.tsv
    done
done
```

```python
import csv, sys

files = ["results_asm5.tsv", "results_asm10.tsv", "results_asm20.tsv", "results_map-ont.tsv"]
readers = [csv.reader(open(f), delimiter='\t') for f in files]
writer = csv.writer(open("results_2026-07.tsv", "wt"), delimiter='\t')

writer.writerow(next(readers[0]))
for r in readers[1:]:
    next(r)

for rows in zip(*readers):
    sums = [sum(map(int, (r[6], r[11], r[12]))) for r in rows]
    idx = min(range(len(sums)), key=lambda i: sums[i])
    writer.writerow(rows[idx])
    print(files[idx])
```

```bash
rm results_asm5.tsv results_asm10.tsv results_asm20.tsv results_map-ont.tsv
```
