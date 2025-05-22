These are the commands I ran for conducting assembler comparison in the Autocycler paper.




# Basecall ONT reads

Did the basecalling on Spartan because it has big GPUs:
```bash
cd /data/scratch/projects/punim1894
mkdir O2024-029; cd O2024-029
scp -r roosta:/home/damg/data/O2024-029/pod5 .
sbatch --job-name=dorado --time=40:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-h100 --gres=gpu:1 --wrap "~/programs/dorado-0.9.5-linux-x64/bin/dorado basecaller --kit-name SQK-RBK114-96 sup pod5 > reads.bam"
```

Demux:
```bash
cd /data/scratch/projects/punim1894/O2024-029
sbatch --job-name=dorado --time=2:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-0.9.5-linux-x64/bin/dorado summary reads.bam > summary.tsv"
sbatch --job-name=dorado --time=2:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-0.9.5-linux-x64/bin/dorado demux --output-dir reads --no-classify -t 8 reads.bam"
```

I'm going to use these five samples from the run, because they are deep, each from a different genus, and (as established in a previous analysis) I know they are uncomplicated (no tricky assembly issues, e.g. heterogeneity):
* AUSMDU00021551, barcode25, _Enterobacter hormaechei_
* AUSMDU00018770, barcode57, _Klebsiella pneumoniae_
* AUSMDU00097349, barcode05, _Listeria welshimeri_
* AUSMDU00026122, barcode63, _Providencia rettgeri_
* AUSMDU00009397, barcode75, _Shigella flexneri_




# Copy reads and basic QC

From this point onward, I did my analysis on the MDU PHL servers.

Create directories:
```bash
mkdir 2025-04_Autocycler_paper
cd 2025-04_Autocycler_paper
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    mkdir "$s"
    mkdir "$s"/reads
    mkdir "$s"/reads_qc
done
```

Copy ONT reads from Spartan and combine into a single file per sample:
```bash
cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei/reads
scp spartan:/data/scratch/projects/punim1894/O2024-029/reads/"*barcode25.bam" .
samtools cat -o nanopore.bam *barcode*.bam

cd ~/2025-04_Autocycler_paper/Klebsiella_pneumoniae/reads
scp spartan:/data/scratch/projects/punim1894/O2024-029/reads/"*barcode57.bam" .
samtools cat -o nanopore.bam *barcode*.bam

cd ~/2025-04_Autocycler_paper/Listeria_welshimeri/reads
scp spartan:/data/scratch/projects/punim1894/O2024-029/reads/"*barcode05.bam" .
samtools cat -o nanopore.bam *barcode*.bam

cd ~/2025-04_Autocycler_paper/Providencia_rettgeri/reads
scp spartan:/data/scratch/projects/punim1894/O2024-029/reads/"*barcode63.bam" .
samtools cat -o nanopore.bam *barcode*.bam

cd ~/2025-04_Autocycler_paper/Shigella_flexneri/reads
scp spartan:/data/scratch/projects/punim1894/O2024-029/reads/"*barcode75.bam" .
samtools cat -o nanopore.bam *barcode*.bam

rm ~/2025-04_Autocycler_paper/*/reads/*barcode*.bam
```

Basic ONT QC (discard reads with a mean qscore <10) and convert to FASTQ:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reads_qc
    samtools fastq -T '*' ../reads/nanopore.bam | tr '\t' ' ' | paste - - - - | awk -F'\t' '{split($1, a, " "); for (i in a) {if (a[i] ~ /^qs:f:/) {split(a[i], b, ":"); if (b[3] >= 10.0) print $0}}}' | sort | tr '\t' '\n' | pigz -p16 > nanopore.fastq.gz
done
```
The above command also sorts the reads by name, which serves to randomise them (good for subsampling later).

Make a 50x-depth ONT read set (mainly for the Unicycler assembly so it doesn't take forever):
```bash
cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei/reads_qc
filtlong --target_bases 269237350 nanopore.fastq.gz | pigz -p16 > nanopore_50x.fastq.gz

cd ~/2025-04_Autocycler_paper/Klebsiella_pneumoniae/reads_qc
filtlong --target_bases 299509800 nanopore.fastq.gz | pigz -p16 > nanopore_50x.fastq.gz

cd ~/2025-04_Autocycler_paper/Listeria_welshimeri/reads_qc
filtlong --target_bases 148627250 nanopore.fastq.gz | pigz -p16 > nanopore_50x.fastq.gz

cd ~/2025-04_Autocycler_paper/Providencia_rettgeri/reads_qc
filtlong --target_bases 223290300 nanopore.fastq.gz | pigz -p16 > nanopore_50x.fastq.gz

cd ~/2025-04_Autocycler_paper/Shigella_flexneri/reads_qc
filtlong --target_bases 241424350 nanopore.fastq.gz | pigz -p16 > nanopore_50x.fastq.gz
```

Copy Illumina reads:
```bash
cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei/reads
cp ~/2024-06_new_assembler_benchmark/AUSMDU00021551/reads/illumina* .

cd ~/2025-04_Autocycler_paper/Klebsiella_pneumoniae/reads
cp ~/2024-06_new_assembler_benchmark/AUSMDU00018770/reads/illumina* .

cd ~/2025-04_Autocycler_paper/Listeria_welshimeri/reads
cp ~/2024-06_new_assembler_benchmark/AUSMDU00097349/reads/illumina* .

cd ~/2025-04_Autocycler_paper/Providencia_rettgeri/reads
cp ~/2024-06_new_assembler_benchmark/AUSMDU00026122/reads/illumina* .

cd ~/2025-04_Autocycler_paper/Shigella_flexneri/reads
cp ~/2024-06_new_assembler_benchmark/AUSMDU00009397/reads/illumina* .
```

Illumina read QC with fastp:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reads_qc
    fastp --in1 ../reads/illumina_1.fastq.gz --in2 ../reads/illumina_2.fastq.gz --out1 illumina_1.fastq.gz --out2 illumina_2.fastq.gz --unpaired1 illumina_u.fastq.gz --unpaired2 illumina_u.fastq.gz
done
```




# Reference genome assembly

To create the reference sequence, I'll use my previously established best practice: Trycycler+Medaka+Polypolish+Pypolca. I'll also do a Unicycler hybrid assembly, which isn't accurate enough for creating the reference sequence but can be useful for checking small plasmids (which can be missed in ONT assemblies).

```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    mkdir reference_assembly
done
```

Unicycler hybrid:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reference_assembly
    if [ -f unicycler_illumina/assembly.gfa ]; then continue; fi
    unicycler -1 ../reads_qc/illumina_1.fastq.gz -2 ../reads_qc/illumina_2.fastq.gz -l ../reads_qc/nanopore_50x.fastq.gz -o unicycler -t 64
done
```

These all completed except for the _Shigella_ genome, where some small plasmids remained tangled together. This was due to my use of Filtlong which left no reads below 15 kbp, so Unicycler couldn't separate those plasmids. But they were easy to separate manually using depth.

Trycycler input assemblies, following the 'extra-thorough' instructions on its wiki:
```bash
declare -A genome_sizes
genome_sizes["Enterobacter_hormaechei"]=5384747
genome_sizes["Klebsiella_pneumoniae"]=5990196
genome_sizes["Listeria_welshimeri"]=2972545
genome_sizes["Providencia_rettgeri"]=4465806
genome_sizes["Shigella_flexneri"]=4828487

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    genome_size=$genome_sizes["$s"]
    
    cd ~/2025-04_Autocycler_paper/"$s"
    mkdir reference_assembly
    cd reference_assembly

    trycycler subsample --reads ../reads_qc/nanopore.fastq.gz --out_dir read_subsets --count 24 --genome_size "$genome_size"
    mkdir assemblies

    for i in 01 07 13 19; do
        canu -p canu -d canu_temp -fast genomeSize="$genome_size" useGrid=false maxThreads=32 -nanopore read_subsets/sample_"$i".fastq
        canu_trim.py canu_temp/canu.contigs.fasta > assemblies/assembly_"$i".fasta
        rm -rf canu_temp
    done

    for i in 02 08 14 20; do
        flye --nano-hq read_subsets/sample_"$i".fastq --threads 32 --out-dir flye_temp
        cp flye_temp/assembly.fasta assemblies/assembly_"$i".fasta
        cp flye_temp/assembly_graph.gfa assemblies/assembly_"$i".gfa
        rm -r flye_temp
    done

    for i in 03 09 15 21; do
        miniasm_and_minipolish.sh read_subsets/sample_"$i".fastq 32 > assemblies/assembly_"$i".gfa
        any2fasta assemblies/assembly_"$i".gfa > assemblies/assembly_"$i".fasta
    done

    for i in 04 10 16 22; do
        ~/programs/NECAT/Linux-amd64/bin/necat.pl config config.txt
        realpath read_subsets/sample_"$i".fastq > read_list.txt
        sed -i "s/PROJECT=/PROJECT=necat/" config.txt
        sed -i "s/ONT_READ_LIST=/ONT_READ_LIST=read_list.txt/" config.txt
        sed -i "s/GENOME_SIZE=/GENOME_SIZE="$genome_size"/" config.txt
        sed -i "s/THREADS=4/THREADS=32/" config.txt
        ~/programs/NECAT/Linux-amd64/bin/necat.pl bridge config.txt
        cp necat/6-bridge_contigs/polished_contigs.fasta assemblies/assembly_"$i".fasta
        rm -r necat config.txt read_list.txt
    done

    for i in 05 11 17 23; do
        echo read_subsets/sample_"$i".fastq > input.fofn
        cp ~/programs/NextDenovo/doc/run.cfg nextdenovo_run.cfg
        sed -i "s/genome_size = 1g/genome_size = "$genome_size"/" nextdenovo_run.cfg
        sed -i "s/parallel_jobs = 20/parallel_jobs = 1/" nextdenovo_run.cfg
        sed -i "s/read_type = clr/read_type = ont/" nextdenovo_run.cfg
        sed -i "s/pa_correction = 3/pa_correction = 1/" nextdenovo_run.cfg
        sed -i "s/correction_options = -p 15/correction_options = -p 32/" nextdenovo_run.cfg
        sed -i "s/-t 8/-t 32/" nextdenovo_run.cfg
        nextDenovo nextdenovo_run.cfg
        cp 01_rundir/03.ctg_graph/nd.asm.fasta nextdenovo_temp.fasta
        rm -r 01_rundir nextdenovo_run.cfg input.fofn
        echo read_subsets/sample_"$i".fastq > lgs.fofn
        cat ~/programs/NextPolish/doc/run.cfg | grep -v "sgs" | grep -v "hifi" > nextpolish_run.cfg
        sed -i "s/parallel_jobs = 6/parallel_jobs = 1/" nextpolish_run.cfg
        sed -i "s/multithread_jobs = 5/multithread_jobs = 32/" nextpolish_run.cfg
        sed -i "s|genome = ./raw.genome.fasta|genome = nextdenovo_temp.fasta|" nextpolish_run.cfg
        sed -i "s|-x map-ont|-x map-ont -t 32|" nextpolish_run.cfg
        nextPolish nextpolish_run.cfg
        cp 01_rundir/genome.nextpolish.fasta assemblies/assembly_"$i".fasta
        rm -r 01_rundir pid*.log.info nextpolish_run.cfg lgs.fofn nextdenovo_temp.fasta
    done

    for i in 06 12 18 24; do
        raven --threads 32 --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_"$i".gfa read_subsets/sample_"$i".fastq > assemblies/assembly_"$i".fasta
    done
done
```

Trycycler cluster:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reference_assembly
    trycycler cluster --assemblies assemblies/*.fasta --reads ../reads_qc/nanopore.fastq.gz --out_dir trycycler --threads 64
done
```

And then Trycycler reconcile:
```bash
trycycler reconcile --reads ../reads_qc/nanopore.fastq.gz --cluster_dir trycycler/cluster_xxx  # run for each good cluster
```

I used the Unicycler assemblies (which are good with small plasmids) to help guide the cluster selection process.

And the remaining Trycycler steps:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reference_assembly
    for c in trycycler/cluster_*; do trycycler msa --threads 64 --cluster_dir "$c"; done
    trycycler partition --reads ../reads_qc/nanopore.fastq.gz --cluster_dirs trycycler/cluster_* --threads 96
    for c in trycycler/cluster_*; do trycycler consensus --cluster_dir "$c"; done
    cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta
done
```

Clean up:
```bash
cd ~/2025-04_Autocycler_paper
rm -r */reference_assembly/read_subsets
rm -r */reference_assembly/assemblies
```

Medaka on each cluster:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reference_assembly
    for c in trycycler/cluster_*; do
        medaka_consensus -i "$c"/4_reads.fastq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r1041_e82_400bps_bacterial_methylation -t 12
        mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
        rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi  # clean up
    done
    cat trycycler/cluster_*/8_medaka.fasta > trycycler_medaka.fasta
    rm trycycler/cluster_*/4_reads.fastq
done
```

Short-read polishing:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reference_assembly
    
    bwa index trycycler_medaka.fasta
    bwa mem -t 24 -a trycycler_medaka.fasta ../reads_qc/illumina_1.fastq.gz > alignments_1.sam
    bwa mem -t 24 -a trycycler_medaka.fasta ../reads_qc/illumina_2.fastq.gz > alignments_2.sam
    polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
    polypolish polish trycycler_medaka.fasta filtered_1.sam filtered_2.sam > trycycler_medaka_polypolish.fasta
    rm *.amb *.ann *.bwt *.pac *.sa *.sam

    pypolca run --careful -a trycycler_medaka_polypolish.fasta -1 ../reads_qc/illumina_1.fastq.gz -2 ../reads_qc/illumina_2.fastq.gz -t 24 -o pypolca
    seqtk seq pypolca/pypolca_corrected.fasta > trycycler_medaka_polypolish_pypolca.fasta
    rm -r pypolca
done
```

Tally up polishing changes:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reference_assembly
    compare_assemblies.py --aligner edlib trycycler.fasta trycycler_medaka.fasta 2> /dev/null | grep -o "\*" | wc -l | tr '\n' '\t'
    compare_assemblies.py --aligner edlib trycycler_medaka.fasta trycycler_medaka_polypolish.fasta 2> /dev/null | grep -o "\*" | wc -l | tr '\n' '\t'
    compare_assemblies.py --aligner edlib trycycler_medaka_polypolish.fasta trycycler_medaka_polypolish_pypolca.fasta 2> /dev/null | grep -o "\*" | wc -l
done
```
Medaka's usefulness was questionable. For the _Enterobacter_ and _Listeria_ genomes, Medaka made one change which Polypolish undid, i.e. the pre-Medaka Trycycler assembly seems to be error-free. And for the _Shigella_ genome, Medaka made multiple changes which Polypolish undid. For the other two (_Klebsiella_ and _Providencia_), Medaka made one change which Polypolish kept.

Dnaapler:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/reference_assembly
    
    dnaapler all -i trycycler_medaka_polypolish_pypolca.fasta -o dnaapler -t 16
    seqtk seq dnaapler/dnaapler_reoriented.fasta > trycycler_medaka_polypolish_pypolca_dnaapler.fasta
    rm -r dnaapler

    cp trycycler_medaka_polypolish_pypolca_dnaapler.fasta ../reference.fasta
done
```

I compared my new assemblies to my old ones, and the only difference was in a very long homopolymer in the _Shigella_ genome. My old assembly had Cx17 but my new one had Cx16. There seem to be no Illumina reads which span this, so I have to lean entirely on ONT reads here. The ONT distribution looks like this:
```
Cx13: 8
Cx14: 20
Cx15: 23
Cx16: 21
Cx17: 24
Cx18: 21
Cx19: 8
Cx20: 6
Cx21: 3
```
So I'll go with the most common length, Cx17, but this one is uncertain.

I then manually renamed the contigs in each `reference.fasta` file to be `chromosome`, `plasmid_1`, `plasmid_2`, etc.

Aside from that long homopolymer in _Shigella_, the reference genome assemblies went smoothly, as I suspected they would (I had previously assembled genomes from this ONT run and selected these five in part because they assembled well). So I'm confident that they are error-free (or close to it).




# Check species

Tool versions used:
* GTDB-Tk 2.4.1
* GTDB R226

To double-check that these genomes match their species labels, I'll run them through GTDB-Tk's classification. Trying with both the Mash-based and classic workflows:
```bash
cd ~/2025-04_Autocycler_paper
mkdir gtdb
mkdir gtdb/in
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_welshimeri Providencia_rettgeri Shigella_flexneri; do
    cp ~/2025-04_Autocycler_paper/"$s"/reference.fasta gtdb/in/"$s".fasta
done
cd gtdb
gtdbtk classify_wf --genome_dir in --out_dir out_1 --cpus 32 --extension fasta --mash_db "$GTDBTK_DATA_PATH"/mash.msh 
gtdbtk classify_wf --genome_dir in --out_dir out_2 --cpus 32 --extension fasta --skip_ani_screen
```

The _Enterobacter_ genome was classified as _Enterobacter hormaechei_C_, i.e. GTDB split _Enterobacter hormaechei_ into multiple species. But that's okay - I'll just leave it as _Enterobacter hormaechei_.

The _Listeria_ genome was a suprise, however: GTDB (both methods) classified it as _Listeria innocua_, not _Listeria welshimeri_. I also checked it against some NCBI genomes to make sure this wasn't a weird GTDB thing, and the results were the same: _Listeria innocua_. So I'll rename this genome:
```bash
cd ~/2025-04_Autocycler_paper
mv Listeria_welshimeri Listeria_innocua
```

Note that the [BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN46906078) for this isolate (which was uploaded previously, not as part of this paper) still says _Listeria welshimeri_.




# Divide ONT reads

The full ONT read set is over 300x deep for each sample, so I'll divide each into 6 non-overlapping 50x deep sets.

```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    mkdir ~/2025-04_Autocycler_paper/"$s"/reads_subsampled
done
```

Run this Python code in each genome's directory, changing the genome size as appropriate:
```python
import gzip
genome_size = 5384747  # Enterobacter_hormaechei
# genome_size = 5990196  # Klebsiella_pneumoniae
# genome_size = 2972545  # Listeria_innocua
# genome_size = 4465806  # Providencia_rettgeri
# genome_size = 4828487  # Shigella_flexneri

target_bp, max_subsets = genome_size * 50, 6
subset, bp = 1, 0
out = open(f"reads_subsampled/{subset}.fastq", "w")
with gzip.open("reads_qc/nanopore.fastq.gz", "rt") as f:
    for rec in zip(*[f]*4):
        _ = out.write("".join(rec))
        bp += len(rec[1].strip())
        if bp >= target_bp:
            if subset < max_subsets:
                out.close()
                subset, bp = subset + 1, 0
                out = open(f"reads_subsampled/{subset}.fastq", "w")
            else:
                break
```

Sanity check - make sure there are no read IDs shared between files:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd  ~/2025-04_Autocycler_paper/"$s"/reads_subsampled
    for f in {1..6}.fastq; do
        cat "$f" | paste - - - - | cut -f1 -d' '
    done | sort | uniq -d
done
```

Sanity check - make each file is 50x depth:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    genome_size=$(seqtk size reference.fasta | cut -f2)
    for f in {1..6}.fastq; do
        reads_size=$(seqtk size reads_subsampled/"$f" | cut -f2)
        awk -v r="$reads_size" -v g="$genome_size" 'BEGIN { printf "%.3fx\n", r/g }'
    done
done
```




# Single-assembler assemblies

Tool versions used:
* Canu 2.3
* Flye 2.9.5
* metaMDBG 1.1
* Miniasm 0.3
* Minipolish 0.1.3
* Necat 0.0.1_update20200803
* NextDenovo 2.5.2
* NextPolish 1.4.1
* Raven-assembler 1.8.3
* Wtdbg 2.5

```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    mkdir ~/2025-04_Autocycler_paper/"$s"/assemblies
done
```

For these assemblies, I'm using my Autocycler helper scripts. These mostly just run the assemblies with default settings, but sometimes do a little bit extra:
* `canu.sh`
  * Runs Canu with `-fast` and otherwise default settings.
  * Then runs my `canu_trim.py` script to trim start-end overlaps and discards 'repeat' and 'bubble' contigs.
* `flye.sh`
  * Runs Flye with `--nano-hq` and default settings.
* `lja.sh`
  * Runs LJA with default settings.
* `metamdbg.sh`
  * Runs metaMDBG with default settings.
  * Then runs my `metamdbg_filter.py` script to remove low-depth contigs. This is needed since metaMDBG is a metagenome assembler and so leaves low-depth contigs in its assemblies.
* `miniasm.sh`
  * Runs miniasm with default settings.
  * Then polishes the result with Minipolish (also default settings).
* `necat.sh`
  * Runs NECAT with default settings.
  * NECAT needs a `config.txt` file for its parameters (not command line arguments), so this script really simplifies running it.
* `nextdenovo.sh`
  * Runs NextDenovo with default settings followed by NextPolish with default settings.
  * NextDenovo and NextPolish need `.fofn` files for their inputs and `.cfg.` files for their parameters (not command line arguments), so this script really simplifies running them.
* `raven.sh`
  * Runs Raven with default settings, but with ` --graphical-fragment-assembly` to additionally save a GFA-format version of the assembly.

```bash
threads=32
jobs=6

declare -A genome_sizes
genome_sizes["Enterobacter_hormaechei"]=5384747
genome_sizes["Klebsiella_pneumoniae"]=5990196
genome_sizes["Listeria_innocua"]=2972545
genome_sizes["Providencia_rettgeri"]=4465806
genome_sizes["Shigella_flexneri"]=4828487

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
        genome_size=$genome_sizes["$s"]
        for assembler in raven redbean miniasm lja metamdbg flye necat nextdenovo canu; do
        for i in {1..6}; do
            echo "nice -n 19 $assembler.sh reads_subsampled/$i.fastq assemblies/${assembler}_$i $threads $genome_size" >> assemblies/jobs.txt
        done
    done
done

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    parallel --jobs "$jobs" --joblog assemblies/joblog.txt --results assemblies/logs < assemblies/jobs.txt
done
```

These were all successful except for one NECAT assembly (`Providencia_rettgeri/assemblies/necat_1.fasta`) which I re-ran with a different thread count:
```bash
cd ~/2025-04_Autocycler_paper/Providencia_rettgeri
necat.sh reads_subsampled/1.fastq assemblies/necat_1 24 4465806
```
This is a common problem I have noticed with NECAT - it just crashes sometimes and using a different number of threads makes it work. I'm actually surprised it only happened once for all 30 NECAT assemblies.




# MAECI assemblies

MAECI is a consensus assembly pipeline: https://github.com/langjidong/MAECI. A bit like Trycycler/Autocycler, but tries to be more end-to-end (e.g. does read QC at the start, polishing at the end).

A bit awkward to install and run. Also, the script seemed to have a bug which I needed to fix in order to make it run: https://github.com/langjidong/MAECI/issues/3

MAECI doesn't have any releases, so I used the current (as of April) pull from GitHub: f1eb3d7. Plus the fix described above.

```bash
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
        mkdir maeci_"$i"
        Multiple-Assembly-Integration.pl -fq reads_subsampled/"$i".fastq -path1 canu -path2 flye -path3 wtdbg2 -path4 racon -genome_size "$genome_size" -outputdir maeci_"$i" -reference reference.fasta -process 32
        cp maeci_"$i"/self-correction/racon.final.fasta assemblies/maeci_"$i".fasta
        rm -rf maeci_"$i"
    done
done
```




# Dragonflye assemblies

[Dragonflye](https://github.com/rpetit3/dragonflye) is bacterial assembly pipeline based on Flye.

I used the latest version as of April 2025: 1.2.1

```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    for i in {1..6}; do
        dragonflye --reads reads_subsampled/"$i".fastq --outdir dragonflye_"$i" --cpus 32
        cp dragonflye_"$i"/contigs.reoriented.fa assemblies/dragonflye_"$i".fasta
        rm -r dragonflye_"$i"
    done
done
```




# Hybracter assemblies

[Hybracter](https://github.com/gbouras13/hybracter) is a bacterial assembly pipeline that includes many steps, including plasmid recovery with Plassembler.

I used the latest version as of April 2025: 0.11.2

```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    for i in {1..6}; do
        hybracter long-single -l reads_subsampled/"$i".fastq -s "$s"_"$i" --auto -o hybracter_"$i" -t 32
        cp hybracter_"$i"/FINAL_OUTPUT/*complete/*_final.fasta assemblies/hybracter_"$i".fasta
        rm -r hybracter_"$i"
    done
done
```

For some reason, Hybracter completely botched the _Shigella_ genome, either producing a terrible fragmented assembly or crashing. I suspect the read QC might be to blame, so I'm trying again with QC turned off:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    for i in {1..6}; do
        hybracter long-single -l reads_subsampled/"$i".fastq -s "$s"_"$i" --auto -o hybracter_"$i" -t 32 --skip_qc
        cp hybracter_"$i"/FINAL_OUTPUT/*complete/*_final.fasta assemblies/hybracter_"$i".fasta
        rm -r hybracter_"$i"
    done
done
```
This worked much better, so I'll use these as my Hybracter assemblies.




# Autocycler assemblies

Used my `autocycler_full.sh` script: https://github.com/rrwick/Autocycler/tree/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick. This includes Plassembler which is given extra clustering weight, especially relevant for the genomes which have small plasmids (_Enterobacter_, _Klebsiella_ and _Shigella_).

I used the latest version as of April 2025: 0.3.1

First, I'll do completely automated assemblies:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for i in {1..6}; do
        cd ~/2025-04_Autocycler_paper/"$s"
        mkdir autocycler_"$i" && cd autocycler_"$i"
        autocycler_full.sh ../reads_subsampled/"$i".fastq 32 4
        cp autocycler_out/consensus_assembly.fasta ../assemblies/autocycler_"$i".fasta
    done
done
```

Make an Autocycler metrics table:
```bash
cd ~/2025-04_Autocycler_paper
autocycler table > metrics.tsv
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for i in {1..6}; do
        autocycler table -a "$s"/autocycler_"$i" -n "$s"_"$i" >> metrics.tsv
    done
done
```

Then I'll make manually curated assemblies. For these, I'll copy the input assemblies I made for the automated Autocycler assemblies:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    for i in {1..6}; do
        mkdir autocycler_"$i"_manual
        mkdir autocycler_"$i"_manual/assemblies
        cp -r autocycler_"$i"/assemblies/*.fasta autocycler_"$i"_manual/assemblies
    done
done
```

And then run the remaining steps for each assembly, doing manual checks along the way:
```bash
autocycler compress -i assemblies -a autocycler_out
autocycler cluster -a autocycler_out

# Manual step: check clustering and re-run with --manual if needed.

for c in autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c"
    if [[ $(wc -c <"$c"/1_untrimmed.gfa) -lt 100000 ]]; then
        autocycler dotplot -i "$c"/1_untrimmed.gfa -o "$c"/1_untrimmed.png
        autocycler dotplot -i "$c"/2_trimmed.gfa -o "$c"/2_trimmed.png
    fi
    autocycler resolve -c "$c"
done
autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

# Manual step: look for non-circularised pieces and re-run earlier steps if needed.
```

Once this are completed, I can then copy final assembly:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    for i in {1..6}; do
        cp autocycler_"$i"_manual/autocycler_out/consensus_assembly.fasta assemblies/autocycler_manual_"$i".fasta
    done
done
```




# Excluded tools

I found some other tools that were vaguely relevant, but they were older and more focused on improving the N50 of a draft assemblies made from short reads (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1213-3). Mainly because of their age (from the 2010s), these programs aren't designed for the same problem as Autocycler (combining a collection of already-finished very good bacterial genome assemblies to a near-perfect bacterial genome assembly).

Some notes on specific tools I tried:
* [quickmerge](https://github.com/mahulchak/quickmerge):
  * Can only merge two assemblies, so would need to be run iteratively in some hierarchical manner to merge more than two.
  * Seems to be focused on improving contiguity, not increasing accuracy, so not really applicable where the input assemblies are already complete.
* [MAC2.0](https://github.com/bioinfomaticsCSU/MAC):
  * Can only merge two assemblies, so would need to be run iteratively in some hierarchical manner to merge more than two.
  * When I ran it on two complete assemblies, it produced an output with two entire copies of the chromosome. So clearly not intended for use on complete bacterial assemblies.
* [Metassembler](https://sourceforge.net/projects/metassembler):
  * Designed for paired-end short read assemblies.




# Count errors and mistakes

Produce the main results table using my `assess_assembly.py` script to check against my ground-truth reference:
```bash
cd ~/2025-04_Autocycler_paper
./assess_assembly.py --header > results.tsv
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/*.fasta; do
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" >> results.tsv
    done
done
```

Also running [Inspector](https://github.com/Maggi-Chen/Inspector) v1.3.1 and [CRAQ](https://github.com/JiaoLaboratory/CRAQ) v1.0.9, both of which assess assemblies using reads (not a ground-truth genome sequence):
```bash
cd ~/2025-04_Autocycler_paper
mkdir inspector_results
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/*.fasta; do
        inspector.py -c "$a" -r "$s"/reads_qc/nanopore.fastq.gz -o inspector -d nanopore -t 32
        name=$(echo $a | sed 's|/assemblies/|_|' | sed 's|.fasta||')
        cp inspector/summary_statistics inspector_results/"$name"
        rm -r inspector
    done
done

cd ~/2025-04_Autocycler_paper
mkdir craq_results
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/*.fasta; do
        craq -g "$a" -sms "$s"/reads_qc/nanopore.fastq.gz -ngs "$s"/reads_qc/illumina_1.fastq.gz,"$s"/reads_qc/illumina_2.fastq.gz --output_dir craq -t 32 -x map-ont
        name=$(echo $a | sed 's|/assemblies/|_|' | sed 's|.fasta||')
        cp craq/runAQI_out/out_final.Report craq_results/"$name"
        rm -r craq
    done
done
```



# Further investigations


## Autocycler sequence-level errors

To allow me to manually inspect the remaining sequence-level errors in Autocycler assemblies, I made Dnaapler-reoriented versions each, and where Dnaapler failed to find a starting position, I manually rotated any sequences to match the reference.
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"/assemblies
    for i in {1..6}; do
        dnaapler all -i autocycler_manual_"$i".fasta -o dnaapler -t 16
        seqtk seq dnaapler/dnaapler_reoriented.fasta > ../autocycler_"$i"_manual/rotated.fasta
        rm -r dnaapler
    done
done
```

I then used my [`compare_assemblies.py`](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki/Comparing-assemblies) script to make a human-readable file of all errors:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    cd ~/2025-04_Autocycler_paper/"$s"
    for i in {1..6}; do
        echo "$s"/autocycler_"$i"_manual/rotated.fasta >> ../autocycler_errors.txt
        compare_assemblies.py autocycler_"$i"_manual/rotated.fasta reference.fasta >> ../autocycler_errors.txt
    done
done
```

Error types, totaled across all 30 Autocycler assemblies:
* Homopolymer-length errors: 67
* Other indels: 13
* Substitutions: 45

Note that these total up to 125 errors, which slightly disagrees with the 123 errors in the results table. This is because the results table was populated by the `assess_assembly.py` script which in a couple cases interpreted a homopolymer-length error near the contig end as missing/extra bases.


## Autocycler structural errors

To check for structural errors in the Autocycler assemblies, I'll run Sniffles to check for structural variants.

Tool versions used:
* Sniffles 2.6.2

```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for i in {1..6}; do
        cd ~/2025-04_Autocycler_paper/"$s"/autocycler_"$i"_manual
        minimap2 -a -x map-ont -t 48 rotated.fasta ../reads_qc/nanopore_50x.fastq.gz | samtools sort > nanopore.bam
        samtools index nanopore.bam
        sniffles -i nanopore.bam -v sniffles.vcf
        rm nanopore.bam*
    done
done
```

All of the Sniffles VCFs were empty, so that's good! But just to check with some positive controls, I'll manually make some structural variants and ensure that Sniffles does find them:
```bash
cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei/autocycler_1_manual
cp rotated.fasta rotated_with_insertion.fasta
cp rotated.fasta rotated_with_deletion.fasta
cp rotated.fasta rotated_with_inversion.fasta

# Manually added structural variants

for sv in insertion deletion inversion; do
    minimap2 -a -x map-ont -t 48 rotated_with_"$sv".fasta ../reads_qc/nanopore_50x.fastq.gz | samtools sort > nanopore.bam
    samtools index nanopore.bam
    sniffles -i nanopore.bam -v sniffles_"$sv".vcf
    rm nanopore.bam*
done
```

Yes, Sniffles found the variant in all three manual cases, so I'm convinced that the negative results from before are genuine.


## Dragonflye error test

I suspect that Dragonflye produced more errors than Flye because of its use of Racon. To test this, I made Dragonflye assemblies of the _Listeria_ genome (smallest so fastest to assemble) with and without Racon:
```bash
cd ~/2025-04_Autocycler_paper/Listeria_innocua
for i in {1..6}; do
    dragonflye --reads reads_subsampled/"$i".fastq --outdir dragonflye_defaults_"$i" --cpus 32
    dragonflye --reads reads_subsampled/"$i".fastq --outdir dragonflye_no_racon_"$i" --cpus 32 --racon 0
done

cd ~/2025-04_Autocycler_paper/Listeria_innocua
../assess_assembly.py --header > dragonflye_results.tsv
for a in dragonflye_*/contigs.reoriented.fa; do
    ../assess_assembly.py -r reference.fasta -a "$a" >> dragonflye_results.tsv
done
```

I was surprised to see that Racon was not to blame! In this test, the no-Racon assemblies were worse than the with-Racon assemblies, so Racon is helping, not hurting. But all assemblies were still much worse than just-Flye assemblies, so there's still something Dragonflye-specific causing problems. Also noteworthy: the error counts in my with-Racon assemblies sometimes differed from my previous runs, so Dragonflye is not deterministic.

New hypothesis: When Dragonflye runs Flye, it uses `-i 0` to turn off Flye polishing. So Dragonflye essentially replaced Flye's internal polishing with Racon-based polishing. I can change this with `--opts '-i 1'` (which actually runs Flye with _both_ `-i 0` and `-i 1` but it seems to work anyway).
```bash
cd ~/2025-04_Autocycler_paper/Listeria_innocua
for i in {1..6}; do
    dragonflye --reads reads_subsampled/"$i".fastq --outdir dragonflye_flye_polish_no_racon_"$i" --cpus 32 --racon 0 --opts '-i 1'
done
for a in dragonflye_flye_polish_no_racon_*/contigs.reoriented.fa; do
    ../assess_assembly.py -r reference.fasta -a "$a" >> dragonflye_results.tsv
done
```

That fixed the errors in 4/6 assemblies, so that seems to be a big part of the problem! But 2/6 assemblies still had problems. Dragonflye's default is to use `--nano-raw` with Flye, but `--nano-hq` is a better choice for modern data. Dragonflye supports this with the `--nanohq` option, so trying that.
```bash
cd ~/2025-04_Autocycler_paper/Listeria_innocua
for i in {1..6}; do
    dragonflye --reads reads_subsampled/"$i".fastq --outdir dragonflye_flye_polish_no_racon_nanohq_"$i" --cpus 32 --racon 0 --opts '-i 1' --nanohq
done
for a in dragonflye_flye_polish_no_racon_nanohq_*/contigs.reoriented.fa; do
    ../assess_assembly.py -r reference.fasta -a "$a" >> dragonflye_results.tsv
done
```

Now the results are pretty much in line with Flye. So I think there are two reasons Dragonflye did worse than Flye:
* It uses Racon to polish instead of Flye's polisher (and it seems as though Flye's polisher is better).
* It uses `--nano-raw` by default but `--nano-hq` does better.

Clean up:
```bash
cd ~/2025-04_Autocycler_paper/Listeria_innocua
rm -r dragonflye_*
```


## Hybracter per-replicon error rates

For the _Enterobacter_ and _Klebsiella_ genomes, Hybracter assemblies often had a noticably higher error rate than Autocycler. I suspect that these errors are coming from the plasmids, because Hybracter uses Plassembler to create its plasmid sequences, and Plassembler relies on Unicycler which can introduce errors into repeat regions. While Medaka could in principle fix these errors, Hybracter seems to only use the Medaka-polished sequence for the chromosome, leaving plasmids in their raw Plassembler form.

To confirm this, I ran the `assess_assembly.py` script on a per-replicon basis for the Hybracter assemblies of the _Enterobacter_ and _Klebsiella_ genomes:
```bash
cd ~/2025-04_Autocycler_paper
./assess_assembly.py --header > hybracter_results.tsv

cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei
grep -A1 "chromosome" reference.fasta > reference_chromosome.fasta
grep -A1 "plasmid_1" reference.fasta > reference_plasmid_1.fasta
grep -A1 "plasmid_2" reference.fasta > reference_plasmid_2.fasta
grep -A1 "plasmid_3" reference.fasta > reference_plasmid_3.fasta
cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei/assemblies
for i in {1..6}; do
    seqtk seq hybracter_"$i".fasta | grep -A1 "chromosome00001" > hybracter_"$i"_chromosome.fasta
    seqtk seq hybracter_"$i".fasta | grep -A1 "plasmid00001" > hybracter_"$i"_plasmid_1.fasta
    seqtk seq hybracter_"$i".fasta | grep -A1 "plasmid00002" > hybracter_"$i"_plasmid_2.fasta
    seqtk seq hybracter_"$i".fasta | grep -A1 "plasmid00003" > hybracter_"$i"_plasmid_3.fasta
    ../../assess_assembly.py  -r ../reference_chromosome.fasta -a hybracter_"$i"_chromosome.fasta >> ../../hybracter_results.tsv
    ../../assess_assembly.py  -r ../reference_plasmid_1.fasta -a hybracter_"$i"_plasmid_1.fasta >> ../../hybracter_results.tsv
    ../../assess_assembly.py  -r ../reference_plasmid_2.fasta -a hybracter_"$i"_plasmid_2.fasta >> ../../hybracter_results.tsv
    ../../assess_assembly.py  -r ../reference_plasmid_3.fasta -a hybracter_"$i"_plasmid_3.fasta >> ../../hybracter_results.tsv
    rm hybracter_"$i"_chromosome.fasta hybracter_"$i"_plasmid_*.fasta
done
cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei
rm reference_chromosome.fasta reference_plasmid_*.fasta

cd ~/2025-04_Autocycler_paper/Klebsiella_pneumoniae
grep -A1 "chromosome" reference.fasta > reference_chromosome.fasta
grep -A1 "plasmid_1" reference.fasta > reference_plasmid_1.fasta
grep -A1 "plasmid_2" reference.fasta > reference_plasmid_2.fasta
grep -A1 "plasmid_3" reference.fasta > reference_plasmid_3.fasta
grep -A1 "plasmid_4" reference.fasta > reference_plasmid_4.fasta
cd ~/2025-04_Autocycler_paper/Klebsiella_pneumoniae/assemblies
for i in {1..6}; do
    seqtk seq hybracter_"$i".fasta | grep -A1 "chromosome00001" > hybracter_"$i"_chromosome.fasta
    seqtk seq hybracter_"$i".fasta | grep -A1 "plasmid00001" > hybracter_"$i"_plasmid_1.fasta
    seqtk seq hybracter_"$i".fasta | grep -A1 "plasmid00002" > hybracter_"$i"_plasmid_2.fasta
    seqtk seq hybracter_"$i".fasta | grep -A1 "plasmid00003" > hybracter_"$i"_plasmid_3.fasta
    seqtk seq hybracter_"$i".fasta | grep -A1 "plasmid00004" > hybracter_"$i"_plasmid_4.fasta
    ../../assess_assembly.py  -r ../reference_chromosome.fasta -a hybracter_"$i"_chromosome.fasta >> ../../hybracter_results.tsv
    ../../assess_assembly.py  -r ../reference_plasmid_1.fasta -a hybracter_"$i"_plasmid_1.fasta >> ../../hybracter_results.tsv
    ../../assess_assembly.py  -r ../reference_plasmid_2.fasta -a hybracter_"$i"_plasmid_2.fasta >> ../../hybracter_results.tsv
    ../../assess_assembly.py  -r ../reference_plasmid_3.fasta -a hybracter_"$i"_plasmid_3.fasta >> ../../hybracter_results.tsv
    ../../assess_assembly.py  -r ../reference_plasmid_4.fasta -a hybracter_"$i"_plasmid_4.fasta >> ../../hybracter_results.tsv
    rm hybracter_"$i"_chromosome.fasta hybracter_"$i"_plasmid_*.fasta
done
cd ~/2025-04_Autocycler_paper/Klebsiella_pneumoniae
rm reference_chromosome.fasta reference_plasmid_*.fasta
```

As expected, the excess errors were often found in large (>100 kbp) plasmids.




# Tarball data for public repo

Assemblies (and references):
```bash
cd ~/2025-04_Autocycler_paper
mkdir assemblies
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    mkdir assemblies/"$s"
    cp "$s"/reference.fasta assemblies/"$s"
    cp "$s"/assemblies/*.fasta assemblies/"$s"
done
find assemblies/ -type f | sort | tar -Jcvf assemblies.tar.xz --owner=0 --group=0 -T -
rm -r assemblies
```

Full reads:
```bash
cd ~/2025-04_Autocycler_paper
mkdir reads_full
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    mkdir reads_full/"$s"
    cp "$s"/reads_qc/illumina_[12].fastq.gz reads_full/"$s"
    cp "$s"/reads_qc/nanopore.fastq.gz reads_full/"$s"
done
find reads_full/ -type f | sort | tar -cvf reads_full.tar --owner=0 --group=0 -T -
rm -r reads_full
```

Subsampled reads:
```bash
cd ~/2025-04_Autocycler_paper
mkdir reads_subsampled
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    mkdir reads_subsampled/"$s"
    cp "$s"/reads_subsampled/*.fastq reads_subsampled/"$s"
done
find reads_subsampled/ -type f | sort | tar -Jcvf reads_subsampled.tar.xz --owner=0 --group=0 -T -
rm -r reads_subsampled
```
