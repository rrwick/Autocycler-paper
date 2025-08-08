These are the commands I ran for conducting all of the analyses in the Autocycler paper.




# Basecall ONT reads

I did the basecalling on [Spartan](https://dashboard.hpc.unimelb.edu.au) because it has big GPUs:
```bash
cd /data/scratch/projects/punim1894
mkdir O2024-029; cd O2024-029
scp -r roosta:/home/damg/data/O2024-029/pod5 .
sbatch --job-name=dorado --time=40:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-h100 --gres=gpu:1 --wrap "~/programs/dorado-0.9.5-linux-x64/bin/dorado basecaller --kit-name SQK-RBK114-96 sup pod5 > reads.bam"
```

Demultiplex the reads:
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
* Flye 2.9.6
* Hifiasm v0.25.0
* LJA v0.2
* metaMDBG 1.1
* Miniasm 0.3, Minipolish 0.1.3
* Myloasm v0.1.0
* Necat 0.0.1_update20200803
* NextDenovo 2.5.2, NextPolish 1.4.1
* Raven 1.8.3
* Wtdbg 2.5

I installed most of these tools with conda. The exception was LJA, where I needed to build from source (see [this GitHub issue](https://github.com/AntonBankevich/LJA/issues/43)).

```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    mkdir ~/2025-04_Autocycler_paper/"$s"/assemblies
done
```

For these assemblies, I'm using [Autocycler helper](https://github.com/rrwick/Autocycler/wiki/Autocycler-helper). It mostly just run the assemblies with default settings, but sometimes does a little bit extra. I also use `--min_depth_rel 0.1` to discard low-depth contigs when that information is available - this is mostly helpful for metaMGBG and Myloasm which are metagenome assemblers and therefore more likely to include low-depth stuff.

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
        for assembler in canu flye hifiasm lja metamdbg miniasm myloasm necat nextdenovo raven redbean; do
        for i in {1..6}; do
            /usr/bin/time -v -o assemblies/"$assembler"_"$i".time autocycler helper "$assembler" --reads reads_subsampled/"$i".fastq --out_prefix assemblies/"$assembler"_"$i" --threads "$threads" --genome_size "$genome_size" --min_depth_rel 0.1
        done
    done
done
```

These were all successful except for two NECAT assemblies (`Enterobacter_hormaechei/assemblies/necat_4.fasta` and `Klebsiella_pneumoniae/assemblies/necat_6.fasta`) which I re-ran with a different thread count:
```bash
cd ~/2025-04_Autocycler_paper/Enterobacter_hormaechei
/usr/bin/time -v -o assemblies/necat_4.time autocycler helper necat --reads reads_subsampled/4.fastq --out_prefix assemblies/necat_4 --threads 24 --genome_size 5384747 --min_depth_rel 0.1

cd ~/2025-04_Autocycler_paper/Klebsiella_pneumoniae
/usr/bin/time -v -o assemblies/necat_6.time autocycler helper necat --reads reads_subsampled/6.fastq --out_prefix assemblies/necat_6 --threads 24 --genome_size 5990196 --min_depth_rel 0.1
```
This is a common problem I have noticed with NECAT - it just crashes sometimes and using a different number of threads makes it work.




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
        /usr/bin/time -v -o assemblies/maeci_"$i".time Multiple-Assembly-Integration.pl -fq reads_subsampled/"$i".fastq -path1 canu -path2 flye -path3 wtdbg2 -path4 racon -genome_size "$genome_size" -outputdir maeci_"$i" -reference reference.fasta -process 32
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
        /usr/bin/time -v -o assemblies/dragonflye_"$i".time dragonflye --reads reads_subsampled/"$i".fastq --outdir dragonflye_"$i" --cpus 32
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
        /usr/bin/time -v -o assemblies/hybracter_"$i".time hybracter long-single -l reads_subsampled/"$i".fastq -s "$s"_"$i" --auto -o hybracter_"$i" -t 32
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
        /usr/bin/time -v -o assemblies/hybracter_"$i".time hybracter long-single -l reads_subsampled/"$i".fastq -s "$s"_"$i" --auto -o hybracter_"$i" -t 32 --skip_qc
        cp hybracter_"$i"/FINAL_OUTPUT/*complete/*_final.fasta assemblies/hybracter_"$i".fasta
        rm -r hybracter_"$i"
    done
done
```
This worked much better, so I'll use these as my Hybracter assemblies.




# Autocycler assemblies

Used my `autocycler_full.sh` script: https://github.com/rrwick/Autocycler/tree/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick. This includes Plassembler which is given extra clustering weight, especially relevant for the genomes which have small plasmids (_Enterobacter_, _Klebsiella_ and _Shigella_).

I used the latest version as of July 2025: 0.5.1

First, I'll do completely automated assemblies. To make the time performance comparable with other tools (which I gave 32 threads), I'm using 8 threads x 4 jobs:
```bash
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for i in {1..6}; do
        cd ~/2025-04_Autocycler_paper/"$s"
        mkdir autocycler_"$i" && cd autocycler_"$i"
        /usr/bin/time -v -o ../assemblies/autocycler_"$i".time autocycler_full.sh ../reads_subsampled/"$i".fastq 8 4
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

Once this is completed, I can then copy the final assembly:
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

I produced the main results table using my `assess_assembly.py` script to check against my ground-truth reference.

Alignment can get a bit fussy when there are indels near the end of an alignment. Also, it is sometimes unclear whether a messy area should be aligned through or clipped off. For these reasons, I ran the script multiple times using different minimap2 preset:
```bash
cd ~/2025-04_Autocycler_paper
./assess_assembly.py --header > results_asm5.tsv
./assess_assembly.py --header > results_asm10.tsv
./assess_assembly.py --header > results_asm20.tsv
./assess_assembly.py --header > results_map-ont.tsv
for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    for a in "$s"/assemblies/*.fasta; do
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset asm5 >> results_asm5.tsv
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset asm10 >> results_asm10.tsv
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset asm20 >> results_asm20.tsv
        ./assess_assembly.py -r "$s"/reference.fasta -a "$a" --alignment_preset map-ont >> results_map-ont.tsv
    done
done
```

I then combined the results tables into one, keeping whichever minimap2 preset produced the fewest total errors (sequence errors + extra bases + missing bases):
```python
import csv, sys

files = ["results_asm5.tsv", "results_asm10.tsv", "results_asm20.tsv", "results_map-ont.tsv"]
readers = [csv.reader(open(f), delimiter='\t') for f in files]
writer = csv.writer(open("results.tsv", "wt"), delimiter='\t')

writer.writerow(next(readers[0]))
for r in readers[1:]:
    next(r)

for rows in zip(*readers):
    sums = [sum(map(int, (r[6], r[11], r[12]))) for r in rows]
    idx = min(range(len(sums)), key=lambda i: sums[i])
    writer.writerow(rows[idx])
    print(files[idx])
```

I also ran [Inspector](https://github.com/Maggi-Chen/Inspector) v1.3.1, which assesses assemblies using long reads, not a ground-truth genome sequence:
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
cd inspector_results
grep "Total small-scale assembly error" *
grep "Structural error" *
```

And [CRAQ](https://github.com/JiaoLaboratory/CRAQ) v1.0.9, which assesses assemblies using reads (both short and long), not a ground-truth genome sequence:
```bash
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
cd craq_results
grep "Genome" * | cut -f6,7 | tr '(' '\t' | tr -d ')'
```

And [BUSCO](https://gitlab.com/ezlab/busco) v6.0.0 (using the closest match lineage for each genome), which assesses assemblies using expected single-copy genes. I used Dnaapler before BUSCO to ensure that genes do not get split across the start-end of circular sequences:
```bash
cd ~/2025-04_Autocycler_paper
mkdir busco_results

declare -A lineage=(
  [Enterobacter_hormaechei]=enterobacter_odb12
  [Klebsiella_pneumoniae]=enterobacteriaceae_odb12
  [Listeria_innocua]=listeria_odb12
  [Providencia_rettgeri]=morganellaceae_odb12
  [Shigella_flexneri]=enterobacteriaceae_odb12
)

for s in Enterobacter_hormaechei Klebsiella_pneumoniae Listeria_innocua Providencia_rettgeri Shigella_flexneri; do
    busco -i "$s"/reference.fasta -m genome --cpu 24 --lineage_dataset ${lineage[$s]} --out_path busco
    cp busco/BUSCO_*/short_summary.*.txt busco_results/"$s"_reference.txt
    cp busco/BUSCO_*/short_summary.*.json busco_results/"$s"_reference.json
    rm -r busco
    for a in "$s"/assemblies/*.fasta; do
        name=$(echo $a | sed 's|/assemblies/|_|' | sed 's|.fasta||')
        dnaapler all -i "$a" -o dnaapler -t 16
        busco -i dnaapler/dnaapler_reoriented.fasta -m genome --cpu 32 --lineage_dataset ${lineage[$s]} --out_path busco
        cp busco/BUSCO_*/short_summary.*.txt busco_results/"$name".txt
        cp busco/BUSCO_*/short_summary.*.json busco_results/"$name".json
        rm -r dnaapler busco
    done
done

cd ~/2025-04_Autocycler_paper/busco_results
for f in $(ls *.txt | grep -v "reference"); do
    awk '/Complete and single-copy BUSCOs/{s=$1} /Complete and duplicated BUSCOs/{d=$1} /Fragmented BUSCOs/{f=$1} /Missing BUSCOs/{m=$1} END{printf "%s\t%s\t%s\t%s\n", s, d, f, m}' "$f"
done
```




# Further investigations

## Autocycler sequence-level errors

To allow me to manually inspect the remaining sequence-level errors in Autocycler assemblies, I made Dnaapler-reoriented versions each, and where Dnaapler failed to find a starting position (some small plasmids), I manually rotated any sequences to match the reference.
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
* Homopolymer-length errors: 76 bp at 57 loci
* Other indels: 15 bp at 15 loci
* Substitutions: 44 bp at 44 loci


## Autocycler structural errors

To check for structural errors in the Autocycler assemblies, I ran Sniffles to check for structural variants.

Tool versions used:
* Sniffles 2.6.3

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

All of the Sniffles VCFs were empty, so that's good! But just to check with some positive controls, I manually made some structural variants to ensure that Sniffles does find them:
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


## Low-depth Autocycler assemblies

To test how well Autocycler handles lower read depths, I used the _Listeria_ genome, for the same reasons I used it in the Dragonflye parameter test: it's small (and therefore fast to assemble) and relatively uncomplicated.

I used the same logic as is in my `autocycler_full.sh` script (used for the main analysis), with these exceptions:
* I used the actual genome size for `autocycler subsample` instead of getting an automatic value with `autocycler helper genome_size`. This is because automatic genome sizes will be unreliable at very low depths.
* With default parameters, Autocycler won't work with read sets shallower than 25x, due to the `--min_read_depth` parameter in `autocycler subsample`. So I used some custom logic here: set `--min_read_depth` to 90% of the read depth but no more than 25x.
* I also used `--max_contigs 1000` with `autocycler compress` and `autocycler cluster` to force it to continue, even when the input assemblies were rubbish. This was to prevent Autocycler from erroring out so I could get a final Autocycler assembly for each read depth.

Run the Autocycler assemblies:
```zsh
mkdir -p ~/2025-04_Autocycler_paper/Listeria_innocua/autocycler_low_depth
cd ~/2025-04_Autocycler_paper/Listeria_innocua/autocycler_low_depth

autocycler table > metrics.tsv

threads=32
jobs=4
read_type=ont_r10

genome_size=2972545
full_read_bases=3347581958
full_read_count=604429

mean_read_length=$(bc -l <<< "$full_read_bases / $full_read_count")
reads_per_1x=$(bc -l <<< "$genome_size / $mean_read_length")

setopt NULL_GLOB

for depth in {01..50}; do
    read_count=$(printf '%.0f' "$(bc -l <<< "$reads_per_1x * $depth")")
    mkdir "$depth"x
    cd "$depth"x
    zcat ../../reads_qc/nanopore.fastq.gz | paste - - - - | shuf | head -n "$read_count" | tr '\t' '\n' > reads.fastq

    typeset -F2 min_read_depth
    (( min_read_depth = depth * 0.9 ))
    (( min_read_depth > 25 )) && (( min_read_depth = 25 ))

    autocycler subsample --reads reads.fastq --out_dir subsampled_reads --genome_size "$genome_size" --min_read_depth "$min_read_depth" 2>> autocycler.stderr
    mkdir -p assemblies
    rm -f assemblies/jobs.txt
    for assembler in raven miniasm flye metamdbg necat nextdenovo plassembler canu; do
        for i in 01 02 03 04; do
            echo "autocycler helper $assembler --reads subsampled_reads/sample_$i.fastq --out_prefix assemblies/${assembler}_$i --threads $threads --genome_size $genome_size --read_type $read_type" --min_depth_rel 0.1 >> assemblies/jobs.txt
        done
    done
    nice -n 19 parallel --jobs "$jobs" --joblog assemblies/joblog.tsv --results assemblies/logs < assemblies/jobs.txt
    for f in assemblies/plassembler*.fasta; do
        sed -i 's/circular=True/circular=True Autocycler_cluster_weight=2/' "$f"
    done
    for f in assemblies/canu*.fasta assemblies/flye*.fasta; do
        sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
    done
    rm subsampled_reads/*.fastq reads.fastq
    autocycler compress -i assemblies -a autocycler_out --max_contigs 1000 2>> autocycler.stderr
    autocycler cluster -a autocycler_out --max_contigs 1000 2>> autocycler.stderr
    for c in autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c "$c" 2>> autocycler.stderr
        autocycler resolve -c "$c" 2>> autocycler.stderr
    done
    autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa 2>> autocycler.stderr
    autocycler table -a . -n "$depth"x >> ../metrics.tsv
    cd ..
done
```

Assess the results:
```bash
cd ~/2025-04_Autocycler_paper/Listeria_innocua/autocycler_low_depth
../../assess_assembly.py --header > results.tsv
for a in */autocycler_out/consensus_assembly.fasta; do
    ../../assess_assembly.py -r ../reference.fasta -a "$a" >> results.tsv
done
```

Looking at the Autocycler logs, I can see that Flye and metaMDBG really carried the weight with lower depth assemblies (e.g. <20x). So I think a low-depth-optimised Autocycler pipeline based on those two assemblers would be possible, and it might (slightly) outperform the generic Autocycler pipeline I used here.


## Mixed/contaminated Autocycler assemblies

For this test, I mixed the _Enterobacter hormaechei_ and _Klebsiella pneumoniae_ genomes at different ratios. I chose these two because they both have high-copy-number small plasmids that will likely appear in the assembly, even at low depths. Also, they are both in the same family (Enterobacteriaceae) which may be close enough to cause issues with some long-read assemblers.

```bash
cd ~/2025-04_Autocycler_paper
mkdir mixed_autocycler_assemblies
cd mixed_autocycler_assemblies
```

This script (`mix_reads.py`) takes care of sampling reads at the appropriate depth:
```python
#!/usr/bin/env python3

import gzip, random, sys

step = int(sys.argv[1])
steps = int(sys.argv[2])

def sample_reads(fastq_gz, read_count):
    with gzip.open(fastq_gz, 'rt') as f:
        reads = [''.join(read) for read in zip(*[f]*4)]
    random.shuffle(reads)
    sys.stdout.write(''.join(reads[:read_count]))

eh_reads = "/home/wickr/2025-04_Autocycler_paper/Enterobacter_hormaechei/reads_qc/nanopore.fastq.gz"
kp_reads = "/home/wickr/2025-04_Autocycler_paper/Klebsiella_pneumoniae/reads_qc/nanopore.fastq.gz"

eh_mean_read_len, kp_mean_read_len = 7048.01591486, 6081.33055366
eh_genome_size, kp_genome_size = 5384747, 5990196

x = step / steps
eh_fraction = x**2 / (x**2 + (1-x)**2)
kp_fraction = 1.0 - eh_fraction
eh_depth = 100.0 * eh_fraction
kp_depth = 100.0 * kp_fraction
eh_read_count = int(round(eh_genome_size * eh_depth / eh_mean_read_len))
kp_read_count = int(round(kp_genome_size * kp_depth / kp_mean_read_len))

sample_reads(eh_reads, eh_read_count)
sample_reads(kp_reads, kp_read_count)
```

Create the read sets:
```bash
for i in {00..50}; do
    cd ~/2025-04_Autocycler_paper/mixed_autocycler_assemblies
    mkdir "$i" && cd "$i"
    ../mix_reads.py "$i" 50 | paste - - - - | sort | tr '\t' '\n' > reads.fastq
done
```

Run Autocycler:
```bash
for i in {00..50}; do
    cd ~/2025-04_Autocycler_paper/mixed_autocycler_assemblies/"$i"
    autocycler_full.sh reads.fastq 32 4
done
```

I assessed the assemblies by categorising each contig as one of the following:
* complete _Enterobacter hormaechei_ sequence
* incomplete _Enterobacter hormaechei_ sequence
* complete _Klebsiella pneumoniae_ sequence
* incomplete _Klebsiella pneumoniae_ sequence

I put this assessment code in `assess_mixed_assembly.py`:
```python
#!/usr/bin/env python3

import sys
import subprocess
import tempfile
from pathlib import Path

def load_fasta(filename, prefix):
    seqs = {}
    with open(filename, "rt") as f:
        while (header := f.readline()):
            seq = f.readline().rstrip("\n")
            name = prefix + header[1:].split()[0]
            seqs[name] = (header.rstrip("\n"), seq)
    return seqs

def write_fasta(filename, seq):
    with open(filename, "wt") as f:
        f.write(f">seq\n{seq}\n")

def mash_identity(seq_1, seq_2):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        fasta_1, fasta_2 = temp_dir / "1.fasta", temp_dir / "2.fasta"
        sketch = temp_dir / "mash.msh"
        write_fasta(fasta_1, seq_1)
        write_fasta(fasta_2, seq_2)
        try:
            subprocess.run(f"mash sketch -k 21 -s 10000 -o {sketch} {fasta_1}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out = subprocess.check_output(f"mash screen {sketch} {fasta_2}", shell=True, text=True, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            return 0.0
        try:
            return float(out.split('\t')[0])
        except ValueError:
            return 0.0

assembly = sys.argv[1]

eh_ref = "/home/wickr/2025-04_Autocycler_paper/Enterobacter_hormaechei/reference.fasta"
kp_ref = "/home/wickr/2025-04_Autocycler_paper/Klebsiella_pneumoniae/reference.fasta"

eh_seqs = load_fasta(eh_ref, "eh_")
kp_seqs = load_fasta(kp_ref, "kp_")
ref_seqs = {**eh_seqs, **kp_seqs}

eh_complete, eh_incomplete, kp_complete, kp_incomplete = 0, 0, 0, 0
for assembly_name, assembly_contig in load_fasta(assembly, "").items():
    assembly_header, assembly_seq = assembly_contig
    best_mash, best_ref, best_complete = 0.0, None, False
    for ref_name, ref_contig in ref_seqs.items():
        ref_header, ref_seq = ref_contig
        complete = ('circular=true' in assembly_header) and (0.9 < (len(assembly_seq) / len(ref_seq)) < 1.1)
        mash = mash_identity(assembly_seq, ref_seq)
        if mash > best_mash:
            best_mash = mash
            best_ref = ref_name
            best_complete = complete
    if best_ref is None:
        continue
    if best_ref.startswith('eh_') and best_complete:
        eh_complete += len(assembly_seq)
    if best_ref.startswith('eh_') and not best_complete:
        eh_incomplete += len(assembly_seq)
    if best_ref.startswith('kp_') and best_complete:
        kp_complete += len(assembly_seq)
    if best_ref.startswith('kp_') and not best_complete:
        kp_incomplete += len(assembly_seq)

print(f"{eh_complete}\t{eh_incomplete}\t{kp_complete}\t{kp_incomplete}")
```

```bash
cd ~/2025-04_Autocycler_paper/mixed_autocycler_assemblies
for i in {00..50}; do
    printf "$i\t"
    ./assess_mixed_assembly.py "$i"/autocycler_out/consensus_assembly.fasta
done
```




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

Low-depth Autocycler assemblies:
```bash
cd ~/2025-04_Autocycler_paper
mkdir low_depth_autocycler
for d in {01..50}; do
    cp ~/2025-04_Autocycler_paper/Listeria_innocua/autocycler_low_depth/"$d"x/autocycler_out/consensus_assembly.fasta low_depth_autocycler/"$d".fasta
done
find low_depth_autocycler/ -type f | sort | tar -Jcvf low_depth_autocycler.tar.xz --owner=0 --group=0 -T -
rm -r low_depth_autocycler
```

Mixed Autocycler assemblies:
```bash
cd ~/2025-04_Autocycler_paper
mkdir mixed_autocycler
for i in {01..50}; do
    cp ~/2025-04_Autocycler_paper/mixed_autocycler_assemblies/"$i"/autocycler_out/consensus_assembly.fasta mixed_autocycler/"$i".fasta
done
find mixed_autocycler/ -type f | sort | tar -Jcvf mixed_autocycler.tar.xz --owner=0 --group=0 -T -
rm -r mixed_autocycler
```
