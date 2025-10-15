# Circos Plot for Wheat Genomic Data

This repository documents how to install and use **Circos** to visualize genomic features (GC content, gene density, custom layers) for wheat genomes on HPC systems like Atlas/CCAST.

---

## 1. Install Circos (Conda Environment)
```bash
ml miniconda3/24.3.0

conda create --prefix /directory/saved/circos_env -c conda-forge circos
source /apps/spack-managed/gcc-11.3.1/miniconda3-24.3.0-tfxfbqlz7yglkzf3fhpokdrwkakbluqw/etc/profile.d/conda.sh
source ~/.bashrc
conda activate /directory/saved/circos_env

or

source /apps/spack-managed/gcc-11.3.1/miniconda3-24.3.0-tfxfbqlz7yglkzf3fhpokdrwkakbluqw/etc/profile.d/conda.sh
conda activate /directory/saved/circos_env
```

## 2. Prepare Configuration Files
### Move into your working directory:
```bash
cd /directory/this/saved/circos_plot/

### Example config (heatmap_circos_wheat.conf):

# ideogram
<ideogram>
<spacing>
default = 0.005r
</spacing>

radius           = 0.90r
thickness        = 20p
fill             = yes
stroke_color     = 127,255,212
stroke_thickness = 2p

show_label       = yes
label_font       = default
label_radius     = dims(image,radius) - 60p
label_size       = 30
label_parallel   = yes
</ideogram>

# ticks
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6
format           = %d

<tick>
spacing        = 10u
size           = 10p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>
</ticks>

karyotype   = x_wheat_karyotype.txt

<image>
dir = /directory/this/saved/circos_plot
file  = circos_4.png
radius         = 1500p
background     = white
angle_offset   = -90
</image>

chromosomes_units = 10000000

<highlights>
<highlight>
file       = x_final_glenn_filtered_output.txt
fill_color = yes
ideogram   = yes
</highlight>
</highlights>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
```

## 3. Run Circos Commands
```bash
circos -conf 1heatmap_circos_wheat.conf
```

## 4. Example Input Data
x_sumai3_gc_content_density_extracted.bed

## 5. Transfer Output from HPC to Local
```bash
scp firstname.lastname@atlas-login.hpc.msstate.edu:/directory/this/saved/circos_plot_sumai3/circos_wheat_17.png C:\Users\firstname.lastname\Documents\Circos_plot\
```

# Data to prepare for Circos layers
To visualize genome assemblies with Circos, all feature data (such as gene density, GC content, repeat content, or custom annotations) must first be extracted and formatted into tab-delimited BED‐like files. Each file should include at minimum the chromosome name, start coordinate, end coordinate, and an optional value or category for coloring or thickness. Separate files should be prepared for each layer you intend to plot (e.g., one file for GC content, another for gene density, another for repeats). These preprocessed files are then referenced in the Circos configuration (.conf) files under <highlights>, <plots>, or <heatmaps> to build multi-layer circular visualizations.

## 1. Karyotype file
The karyotype file lists each chromosome (ta1A–ta7D) with its size and an RGB color. This file is referenced in your Circos configuration under `karyotype = x1-wheat-karyotype`.
```bash
- **Columns:**  
  - `chr` – keyword indicating a chromosome line  
  - `-` – placeholder  
  - `ta1A` – internal ID used in Circos  
  - `1A` – label displayed in the plot  
  - `0` – start coordinate  
  - `603896264` – chromosome length (bp)  
  - `137,236,218` – RGB color  

This file must be saved in your Circos plot directory, e.g.:  
`/directory/this/saved/circos_plot_wheat/x1-wheat-karyotype`.
```

## 2. Subtelomeric Tandem Repeats
This workflow searches for the **A6-10 subtelomeric tandem repeat** (AY249980.1) in wheat assemblies using **BLASTN**, then prepares **Circos highlights** files (BED-like with `fill_color=`) for plotting.

### Setup & Query FASTA
```bash
# Work dir
mkdir -p /directory/this/saved/subtelomeric-tandem-repeats
cd /directory/this/saved/subtelomeric-tandem-repeats

# Create query FASTA (AY249980.1 Aegilops tauschii clone A6-10)
nano AY249980.1.fasta

>AY249980.1 Aegilops tauschii clone A6-10 subtelomeric tandem repeat sequence
AGGTCGACGGTATCCATAAGCTTCCAAAGGATAATGAATTGCCCGACAGTACGTACGCAGCAAAGAAGGT
CGTTTGCCCACTATGATTGGAGGTGCAGAACATACATGCATGCCCTAATGACTGCATCCTCTACCGCGAT
GCGTACAAGGATTTGAACGCATGCCCGGTATGCGGTGCATTGCGGTATAAGATCAGACGAGGTGACCCTG
GTGATGTTAATGGCGAGCCCCGCCAGGAAGAGGGTTCCTGCGAAGGTGATGTGGTATGCTCCTATAATAC
CACGGCTGAAACGACTGTTCACAAACAAAAAGCATGCCAAGTTGATGCGATGGCACCGTGAGGACCGTAA
GAAAGACGGGAAGTTGAGAGCACCCGCTGACGGGTCATAGTGTACAAAAATCTAGAGAAAGTACTGGGCT
GAGTTTGCAGGTGACCCAAGGAACATATGGTTTGCTTTAAGCGCGGATGGCATTAATCCTTTCGGGGAGC
AGAGCAGCGATCACAACACTTGGCCCGTGACTCTATGTATGTATAACCTTCCTCCTTGGATGTGCATGAA
GCGGAAGTTCATTATGATGACAATTTTCATCCAAGGCCCTAAGCAACCCGGCAACGACATTGATGTGTAC
CTAAGGGCATTAGTTGAAGAACTTCTACAGTTGCGGAATGGAAACCGTGTACGTGCATGGGATGAGCACA
GACTGGAGGAATTTTACCTTGCACGATTGCTGTTTGTAACATCCATGATTGGCCCGCTCTCATAACCCTT
CCGGACACACAACAAGGGGTGCCCCGCATGCACCCACTGCTTACTTG
```

### BLAST+: Make Databases (per genome) & Search
```bash
ml blast-plus/2.14.1

# Make BLAST databases
makeblastdb -in /where/this/saved/Wheat_pm_v1.fasta  -dbtype nucl -out Wheat_db

# Run BLASTN (tabular outfmt 6)
blastn -query AY249980.1.fasta -db Wheat_db  -out wheat_A6-10_hits.txt  -evalue 1e-10 -outfmt 6
```

### Filter High-Identity Hits at Chromosome Ends
```bash
# Filter 91–100% identity AND within first/last 100 Mb (Wheat chromosome sizes used below)
awk '$3 >= 91 && $3 <= 100 && \
(($2 == "chr1A" && ($9 <= 100000000 || $9 >= 500972114)) || \
 ($2 == "chr1B" && ($9 <= 100000000 || $9 >= 612208373)) || \
 ($2 == "chr1D" && ($9 <= 100000000 || $9 >= 408628255)) || \
 ($2 == "chr2A" && ($9 <= 100000000 || $9 >= 693503254)) || \
 ($2 == "chr2B" && ($9 <= 100000000 || $9 >= 707260945)) || \
 ($2 == "chr2D" && ($9 <= 100000000 || $9 >= 564395422)) || \
 ($2 == "chr3A" && ($9 <= 100000000 || $9 >= 661249593)) || \
 ($2 == "chr3B" && ($9 <= 100000000 || $9 >= 770350482)) || \
 ($2 == "chr3D" && ($9 <= 100000000 || $9 >= 540418347)) || \
 ($2 == "chr4A" && ($9 <= 100000000 || $9 >= 663779771)) || \
 ($2 == "chr4B" && ($9 <= 100000000 || $9 >= 601929352)) || \
 ($2 == "chr4D" && ($9 <= 100000000 || $9 >= 435113464)) || \
 ($2 == "chr5A" && ($9 <= 100000000 || $9 >= 619264322)) || \
 ($2 == "chr5B" && ($9 <= 100000000 || $9 >= 643323281)) || \
 ($2 == "chr5D" && ($9 <= 100000000 || $9 >= 485869797)) || \
 ($2 == "chr6A" && ($9 <= 100000000 || $9 >= 527955843)) || \
 ($2 == "chr6B" && ($9 <= 100000000 || $9 >= 642456823)) || \
 ($2 == "chr6D" && ($9 <= 100000000 || $9 >= 408367107)) || \
 ($2 == "chr7A" && ($9 <= 100000000 || $9 >= 649415582)) || \
 ($2 == "chr7B" && ($9 <= 100000000 || $9 >= 668511397)) || \
 ($2 == "chr7D" && ($9 <= 100000000 || $9 >= 558343710))) \
{print $2, $9, $10}' wheat_A6-10_hits.txt > wheat_filtered_output_100Mb.txt
```

### Rename chr* → ta* and Add Circos fill_color=
```bash
# Rename chromosomes to match Circos karyotype IDs (chr* -> ta*)
sed -e 's/chr1A/ta1A/g' -e 's/chr1B/ta1B/g' -e 's/chr1D/ta1D/g' \
    -e 's/chr2A/ta2A/g' -e 's/chr2B/ta2B/g' -e 's/chr2D/ta2D/g' \
    -e 's/chr3A/ta3A/g' -e 's/chr3B/ta3B/g' -e 's/chr3D/ta3D/g' \
    -e 's/chr4A/ta4A/g' -e 's/chr4B/ta4B/g' -e 's/chr4D/ta4D/g' \
    -e 's/chr5A/ta5A/g' -e 's/chr5B/ta5B/g' -e 's/chr5D/ta5D/g' \
    -e 's/chr6A/ta6A/g' -e 's/chr6B/ta6B/g' -e 's/chr6D/ta6D/g' \
    -e 's/chr7A/ta7A/g' -e 's/chr7B/ta7B/g' -e 's/chr7D/ta7D/g' \
    wheat_filtered_output_100Mb.txt > wheat_filtered_output_100Mb.txt

# Add Circos highlight color (teal = 0,128,128)
awk '{print $1, $2, $3, "fill_color=0,128,128"}' updated-wheat-filtered-output-100Mb.txt \
  > updated-wheat-filtered-output-100Mb.txt

# (Optional) Simpler (no end-zone filter) path:
# sed ... wheat_filtered_output.txt > updated_wheat_filtered_output.txt
# awk '{print $1, $2, $3, "fill_color=0,128,128"}' updated-wheat-filtered-output-100Mb.txt > x2-subtelo-wheat-filtered-output-100Mb
```

### Move Highlights File to Circos Project
```bash
cp x2-subtelo-wheat-filtered-output-100Mb /directory/this/saved/circos_plot_wheat/
```
In your Circos .conf:
```bash
<highlights>
<highlight>
file       = x2-final-wheat-filtered-output-100Mb
fill_color = yes
ideogram   = yes
</highlight>
</highlights>
```

## 3. Centromere-Like Tandem Repeats (Aegilops tauschii 6C6-3/6C6-4) — BLASTN & Circos Prep

This workflow finds **centromere-associated tandem repeats** in wheat assemblies using **BLASTN** and prepares **Circos** highlight files (BED-like with `fill_color=`).

---

### Combine Centromeric Repeat Sequences
```bash
mkdir -p /directory/this/saved/centromeric-tandem-repeats
cd /directory/this/saved/centromeric-tandem-repeats

nano combined-fasta-centromere

>AY249981.1 Aegilops tauschii clone 6C6-4 centromere-specific tandem repeat sequence
TTCGCCCTTGAATGGACAACAGCATTTAAAACTAAATTTGGTTTATATGAGTGGTTAGTCATGCCTTTTG
GGTTAACTAATGCACCTAGTACTTTCATGAGACTAATGAACGAAGTTTTACGTGCTTTCATTGGATGATT
TGTGGTAGTCTATTGTGATGATATACTGATTTATAGCAAATCTTTGGAAGAACATTTGGAACATTTACGC
GCTGTTTTTATTGCTCTACGTGATGCATGTTTGTTTGGTAACCTTGGGAAGTGCACCTTTTGCACCGACC
GAGTATCTTTTCTTGGCTATGTTGTTACTCCACAGGGAATTCAAGTTGGTAAAGCCAAGATTGAAGCTAT
TGAGAGTTGGCCGCAGCCCAAAACGGTCACACAAGTGCGGAGTTTTCTTGGACTCGCTGGTTTCTATAGG
CGTTTTGTGAGAGATTTTAGCACCATCGCTGCACCTCTCAACGAGCTTACAAAGAAGGATGTGCCTTTTG
TTTGGGGTACCGCACAGGAAGAAGCCTTCACGGTATTGAAAGATAAGTTGACACATGCTCCTTTACTCCA
ACTTCCTGATTTTAATAAGACTTTTGAGCTTGAATGTGATGCTAGTGGCATTGGATTAGGAGGTGTGTTA
TTACAAGATGGCAAACCTGTTGCATATTTTTCTGAAAAATTGAGTGGGCCTAGTCTGAATTATTCTACTT
ATGATAAAGAACTATATGCTCTTGTTCGGACTTTAGAAACATGGCAACATTATTTATGGCCCAAAGAATT
TGTTATACATTCTGATCATGAATCTTTGAAACATATTAAAAGTCAAGCTAAACTGAATCGTAGACATGCT
AAATGGGTTGAATTCATTGAAACTTTCCCTTATGTCATTAAACACAAGAAGGGTAAAGAAAATGTTATTG
CTGATGCATTGTCTCGTCGCTATACTATGCTTTCACAAATTGACTTCAAAATATTTGGTTTGGAGACCAT
TAAAGATCAATATGTGCATGATGCTGATTTTAAAGATGTATTGCAGAATTGTAAAGAAGGATGAACCTGG
AAACAGTTTGTCG
>AY249982.1 Aegilops tauschii clone 6C6-3 centromere-specific tandem repeat sequence
ATTCGCCCTTACCCGCCTGCATCTTCTACTTCCACTGCACCAGCACATCCATCAGGTTCCACCTCTAGCC
GTGATACAAGAAAGCAGGCACAACCACCACTATCTGCCAAGAGCGCACCTGCCGGGCCTGCACAGAGCTC
TTCTTCTTCCATGGCACCGACAGGGCACACAAGTGATATTATTTGTCGTCGTTGTAAGGGAAGAGGACAT
TATGCGAGAGAATGCAAATCTCAGCGTGTAATGATTGCTACTGAGGATGGTGGGTATGAGTCCGCTAGTG
ACTATGATGAGGAGACTTTGGCTCTTATTACACGTGAAGAACATGGTGGAGATGATTCTGAAAATGAGAC
GCAATACATGGCTCCTGAAGACGCTGACAGGTATGAATGTTTAGTTGCTCAACGTGTTTTGAGTGTGCAG
GTCACACAAGCAGAGCAAAATCAGAGGCATAATTTGTTCCATACAAAGGGAGTTGTGAAGGAACGTTCTG
TTCGCGTCATCATAGATGGAGGGAGCTGTAACAACTTGGCTAGCATGGAGATGGTGGAGAAGCTATCTCT
CACCACAAGACCACATCCACATCCTTACTACATCCAACGGTTCAACAACAGCGGCAAGGTTAAGGTAACA
CGTACTGTTCGTGTGCATTTTAGTATCTCTACATATGCTGATTATGTTGATTGTGATGTGGTACCCATGC
AAGCATGTTCCTTATTACTTGGTAGACCATGGCAATTTGATAAAAATTCTGTACACCATGGTAGAAACAA
TCAGTATACTCTTGTTCATAAGGATAAAAATATTACTTTGCTTCCTATGACTCCTGATTCCATTTTGAAA
GATGACATTAATAGAGCTAATAAAGCAAAAACAGGAGAAAAATAAGAGTGAAAATCAGATTGTGGCAAAA
GAATTTGAGCAACAAATGAAACTTAATAATAAACCATCTAGTGTTGTTTCTGAAATTAAATTGAAAAGTG
CATGTTTACTTGCCACAAAATCTGATATTGATGAGCTAGATTTCAGCAAATCTGTTTGCTATGCTTTTGT
GTGCAAAGTAGGGCGATTATTTTAAGGGGCGGAA
```

### BLAST+: Search Against Each Genome
Reuse the databases you already made (Wheat_db). If needed, create them with makeblastdb first.
```bash
ml blast-plus/2.14.1

cd /directory/this/saved/centromeric-tandem-repeats

blastn -query combined-fasta-centromere -db Wheat_db  -out combined-centromere-wheat_hits.txt  -evalue 1e-10 -outfmt 6 -num_threads 20
```

### Filter Specific Hits on Wheat chromosome
### Extract % Identity Hits
```bash
# Keep exact matches (100% identity); output: query_id, subject_chr, sstart, send, evalue
awk '$3 == 100 {print $1, $2, $9, $10, $11}' combined-centromere-wheat_hits.txt > centromere-wheat-100percent-hits

# For 95–100% identity
awk '$3 >= 95 && $3 <= 100 {print $1, $2, $9, $10, $11}' combined-centromere-wheat_hits.txt > centromere-wheat-95-100percent-hits.txt
```

### Print lines for chromosomes 
```bash
awk '$2 == "chr1A"' centromere-wheat-95-100percent-hits.txt
```

### Counts how many coordinates start with each prefix e.g. 20–29
```bash
awk '
$2=="chr1A" && $3>=95 && $3<=100 && $4>1000 {
  for (i=9; i<=10; i++) {
    if ($i ~ /^2[0-9]/) {
      prefix = substr($i,1,2)
      count[prefix]++
    }
  }
}
END {
  print "---- Counts by prefix ----"
  for (p=20; p<=29; p++) {
    printf "%2d*: %d\n", p, count[p]+0
  }
}' combined-centromere-wheat_hits.txt
```

### To identify only the highest and lowest positions for Wheat on chrromosome with high identity (95–100%), length >1000, and coordinates beginning with n* or n* Mb:
```bash
awk '
$2=="chr1A" && $3>=95 && $3<=100 && $4>1000 {
  if ($9 ~ /^21[0-9]/ || $10 ~ /^21[0-9]/) {
    if (min=="" || $9<min) min=$9
    if (min=="" || $10<min) min=$10
  }
  if ($9 ~ /^22[0-9]/ || $10 ~ /^22[0-9]/) {
    if (max=="" || $9>max) max=$9
    if (max=="" || $10>max) max=$10
  }
}
END {print "Lowest (21*):", min; print "Highest (22*):", max}' combined-centromere-wheat_hits.txt
```

### Rename chr* → ta* for Circos Karyotype IDs
```bash
sed -e 's/chr1A/ta1A/g' -e 's/chr1B/ta1B/g' -e 's/chr1D/ta1D/g' \
    -e 's/chr2A/ta2A/g' -e 's/chr2B/ta2B/g' -e 's/chr2D/ta2D/g' \
    -e 's/chr3A/ta3A/g' -e 's/chr3B/ta3B/g' -e 's/chr3D/ta3D/g' \
    -e 's/chr4A/ta4A/g' -e 's/chr4B/ta4B/g' -e 's/chr4D/ta4D/g' \
    -e 's/chr5A/ta5A/g' -e 's/chr5B/ta5B/g' -e 's/chr5D/ta5D/g' \
    -e 's/chr6A/ta6A/g' -e 's/chr6B/ta6B/g' -e 's/chr6D/ta6D/g' \
    -e 's/chr7A/ta7A/g' -e 's/chr7B/ta7B/g' -e 's/chr7D/ta7D/g' \
    centromere-wheat-100percent-hits > updated-centromere-wheat-100percent-hits
```

### Prepare your data for layer 2
```bash
# Example of the file
nano x3-centromere-glenn-highlights

ta1A	0	201260369	fill_color=aquamarine
ta1A	201260370	219476881	fill_color=teal
ta1A	219476882	603896264	fill_color=aquamarine
ta1B	0	254692168	fill_color=aquamarine
ta1B	254692169	264373010	fill_color=teal
ta1B	264373011	710210752	fill_color=aquamarine
ta1D	0	230211361	fill_color=aquamarine
```

### Keep Only chr/ta & Coordinates; Add Circos Color
```bash
# Keep subject chr (now ta*), start, end
awk '{print $2, $3, $4}' updated-centromere-wheat-100percent-hits \
  > x-centromere-wheat-100percent-hits

# Add a color tag for Circos highlights (e.g., orange)
awk '{print $1, $2, $3, "fill_color=176,72,103"}' x-centromere-wheat-100percent-hits \
  > x3-centromere-wheat-highlights
```

### Move Highlights File to Circos Project
```bash
cp x3-centromere-wheat-highlights /directory/this/saved/circos_plot_wheat/
```
In your Circos .conf:
```bash
<highlights>
<highlight>
file       = x3-centromere-wheat-highlights
fill_color = yes
ideogram   = yes
</highlight>
</highlights>
```

## 4. Noncoding gene dentsity (infernal RNAs)
### Assuming you already generated updated_wheat_ncrna.bed (chrom  start  end),
### this command aggregates non-coding RNA counts into 100 kb windows for Circos:

```bash
awk -v bin_size=100000 '
{
  start_bin = int($2 / bin_size) * bin_size
  end_bin   = int($3 / bin_size) * bin_size
  for (bin = start_bin; bin <= end_bin; bin += bin_size) {
    density[$1][bin]++
  }
}
END {
  for (chrom in density) {
    for (bin in density[chrom]) {
      print chrom, bin, bin + bin_size - 1, density[chrom][bin]
    }
  }
}' updated_wheat_ncrna.bed | sort -k1,1 -k2,2n > x4-ncrna-glenn-100kb-windows

# Example file
ta1A 0 99999 33
ta1A 100000 199999 34
ta1A 200000 299999 6
```

## 5. Noncoding gene density (infernal RNAs)
### Assuming you already generated wheat.fasta.2.7.7.80.10.50.500.dat,
### this command prepare BED Files for Genome-wide Density (Circos):

```bash
awk '/Sequence:/ {chr=$2} /^[0-9]/ {print chr"\t"($1-1)"\t"$2}' wheat.fasta.2.7.7.80.10.50.500.dat > wheat_tandem_repeats.bed

# Make windows:
ml bedtools2/2.31.1
bedtools makewindows -g wheat.fasta.fai -w 1000000 > wheat_genome_1Mb_windows.bed

# Coverage:
bedtools coverage -a wheat_genome_1Mb_windows.bed -b wheat_tandem_repeats.bed > wheat_tandem_repeat_density.txt

# Format for Circos:
awk '{gsub(/^chr/, "ta", $1); print $1, $2, $3, $7}' wheat_tandem_repeat_density.txt > x5_tandem_repeat_density_wheat
```

# 6. LTR density
### Assuming you already generated wheat.fasta.finder.combine.gff3 from LTR_FINDER_parallel,
### this command prepare 1 Mb window density (bp per bin) for long_terminal_repeat features:

```bash
awk -v window=1000000 '
$3 == "long_terminal_repeat" {
    chr = $1
    start = $4
    end = $5
    bin = int(start / window) * window
    density[chr, bin] += (end - start + 1)
}
END {
    for (key in density) {
        split(key, keys, SUBSEP)
        chr = keys[1]
        bin_start = keys[2]
        bin_end = bin_start + window
        print chr, bin_start, bin_end, density[key]
    }
}' ltr-wheat.gff3 > ltr-wheat-density-circos.txt

# Rename chr* to ta* (e.g., chr1A -> ta1A)
sed -E 's/chr([1-7][ABD])/ta\1/' ltr-wheat-density-circos.txt > x6-ltr-rollag-density

# Get min/max density
awk 'BEGIN {max=0; min=1e18} {if ($4>max) max=$4; if ($4<min) min=$4} END {print "Max Density:", max; print "Min Density:", min}' x6-ltr-rollag-density
```

Maintainer:

Ruby Mijan



