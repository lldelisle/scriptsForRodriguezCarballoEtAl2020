# scriptsForRodriguezCarballoEtAl2020
all the scripts needed to reproduce the figures in Rodriguez-Carballo et al. 2020 from raw data.
 

The full pipeline is described in the [pipeline.sh](https://github.com/lldelisle/scriptsForRodriguezCarballoEtAl2020/blob/master/pipeline.sh) file.

The slurm_scripts directory contains all scripts which are launched during the pipeline on the slurm job scheduler (for the publication it was launched on [SCITAS](http://scitas.epfl.ch/)).

The script directory contains all R, python, or bash scripts used in the slurm scripts.

The tables directory contains all tables, annotations, text files used in the slurm scripts.

The outputs directory contains all annotation files, ini files needed to reproduce the figures as well as the outputs of pyGenomeTracks which were then modified in illustrator to make the final figures.
