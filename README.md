# Circos Plot for Wheat Genomic Data

This repository documents how to install and use **Circos** to visualize genomic features (GC content, gene density, custom layers) for wheat genomes on HPC systems like Atlas/CCAST.

---

## 1. Install Circos (Conda Environment)
```bash
conda create --prefix /directory/saved/circos_env -c conda-forge circos
source ~/.bashrc
conda activate /directory/saved/circos_env
```

## 2. Prepare Configuration Files
### Move into your working directory:
```bash
cd /directory/this/saved/circos_plot/

### Example config (6heatmap_circos_wheat.conf):

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

Maintainer:
Ruby Mijan



