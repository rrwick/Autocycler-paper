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
        /usr/bin/time -v -o assemblies/myloasm_v0.2.0_"$i".time autocycler helper myloasm --reads reads_subsampled/"$i".fastq --out_prefix assemblies/myloasm_v0.2.0_"$i" --threads "$threads" --genome_size "$genome_size" --min_depth_rel 0.1
        /usr/bin/time -v -o assemblies/metamdbg_v1.2_"$i".time autocycler helper metamdbg --reads reads_subsampled/"$i".fastq --out_prefix assemblies/metamdbg_v1.2_"$i" --threads "$threads" --genome_size "$genome_size" --min_depth_rel 0.1
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
    for a in "$s"/assemblies/*_v*.fasta; do
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
writer = csv.writer(open("results_new_versions.tsv", "wt"), delimiter='\t')

writer.writerow(next(readers[0]))
for r in readers[1:]:
    next(r)

for rows in zip(*readers):
    sums = [sum(map(int, (r[6], r[11], r[12]))) for r in rows]
    idx = min(range(len(sums)), key=lambda i: sums[i])
    writer.writerow(rows[idx])
    print(files[idx])
```

[Inspector](https://github.com/Maggi-Chen/Inspector) v1.3.1:
```bash
cd ~/2025-04_Autocycler_paper
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/*_v*.fasta; do
        inspector.py -c "$a" -r "$s"/reads_qc/nanopore.fastq.gz -o inspector -d nanopore -t 32
        name=$(echo $a | sed 's|/assemblies/|_|' | sed 's|.fasta||')
        cp inspector/summary_statistics inspector_results/"$name"
        rm -r inspector
    done
done
cd ~/2025-04_Autocycler_paper/inspector_results
grep "Total small-scale assembly error" *_v*
grep "Structural error" *_v*
```

[CRAQ](https://github.com/JiaoLaboratory/CRAQ) v1.0.9:
```bash
cd ~/2025-04_Autocycler_paper
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/*_v*.fasta; do
        craq -g "$a" -sms "$s"/reads_qc/nanopore.fastq.gz -ngs "$s"/reads_qc/illumina_1.fastq.gz,"$s"/reads_qc/illumina_2.fastq.gz --output_dir craq -t 32 -x map-ont
        name=$(echo $a | sed 's|/assemblies/|_|' | sed 's|.fasta||')
        cp craq/runAQI_out/out_final.Report craq_results/"$name"
        rm -r craq
    done
done
cd ~/2025-04_Autocycler_paper/craq_results
grep "Genome" *_v* | cut -f6,7 | tr '(' '\t' | tr -d ')'
```

[BUSCO](https://gitlab.com/ezlab/busco) v6.0.0:
```bash
cd ~/2025-04_Autocycler_paper

declare -A lineage=(
  [Enterobacter_hormaechei]=enterobacter_odb12
  [Klebsiella_pneumoniae]=enterobacteriaceae_odb12
  [Listeria_innocua]=listeria_odb12
  [Providencia_rettgeri]=morganellaceae_odb12
  [Shigella_flexneri]=enterobacteriaceae_odb12
)

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/*_v*.fasta; do
        name=$(echo $a | sed 's|/assemblies/|_|' | sed 's|.fasta||')
        dnaapler all -i "$a" -o dnaapler -t 16
        busco -i dnaapler/dnaapler_reoriented.fasta -m genome --cpu 32 --lineage_dataset ${lineage[$s]} --out_path busco
        cp busco/BUSCO_*/short_summary.*.txt busco_results/"$name".txt
        cp busco/BUSCO_*/short_summary.*.json busco_results/"$name".json
        rm -r dnaapler busco
    done
done

cd ~/2025-04_Autocycler_paper/busco_results
for f in $(ls *_v*.txt); do
    printf "$f\t"
    awk '/Complete and single-copy BUSCOs/{s=$1} /Complete and duplicated BUSCOs/{d=$1} /Fragmented BUSCOs/{f=$1} /Missing BUSCOs/{m=$1} END{printf "%s\t%s\t%s\t%s\n", s, d, f, m}' "$f"
done
```
