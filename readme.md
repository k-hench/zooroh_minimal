# Minimal example of a ZooRoH run

produces two output files:

```
results/zooroh/
├── inbreeding_by_HBD_class.tsv
└── roh_segments.tsv.gz
```

## setup

Requires `snakemake`, the bash variable `${CDATA}` to be set,  one container ([r_zooroh_v0.2](https://hub.docker.com/r/khench/r_zooroh/tags), located at `${CDATA}/apptainer_local`), one `conda` environment (`popgen_basics`).

pull container:

```sh
cd $CDATA/apptainer_local
apptinaer pull docker://khench/r_zooroh:v0.2
```

install `conda` env (if needed)

```sh
cd code/workflow/envs
mamba env create -f popgen_basics.yml
```

plug in data into the `data` folder (search at the elephant_seal stash...)

## run

execute the pipeline from the `code` directory

```sh
cd code
snakemake -c 7 --use-conda --use-singularity --singularity-args "--bind $CDATA"
```