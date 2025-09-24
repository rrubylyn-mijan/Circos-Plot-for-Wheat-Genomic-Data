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


Maintainer:
Ruby Mijan



