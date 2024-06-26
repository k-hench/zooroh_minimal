#!/usr/bin/env Rscript
# run:
# Rscript <gen_file> <sample_file> <out_dir> <out_roh> <out_f_tab> <gt_fmt> <rate_type> <n_threads>
# eg:
# Rscript \
#   results/genotypes/mirang_test.gen.gz \
#   results/mirang_test.samples \
#   results/zooroh \
#   results/zooroh/roh_segments.tsv.gz \ 
#   results/zooroh/inbreeding_by_HBD_class.tsv \
#   gp bp 7
# ---------------------------------
# snakemake ins-and-outs:
# input:
#   beds = expand( "../results/recombination/fasteprr_4/pop_recrate_{mscaf}.bed.gz", mscaf = AUTOSOMES )
# output:
#   pdf_r = "../results/img/.pdf",
#   r_plot_r = "../results/img/R/.Rds"
# ---------------------------------
# interactive quick-start copy-and-paste template
# args <- c("results/genotypes/mirang_test.gen.gz",
#           "results/mirang_test.samples",
#           "results/zooroh",
#           "results/zooroh/roh_segments.tsv.gz",
#           "results/zooroh/inbreeding_by_HBD_class.tsv",
#            "gp", "bp", 7)
# ---------------------------------

# --- overhead / logistics ---
library(tidyverse)
library(here)
library(glue)
library(RZooRoH)

args <- commandArgs(trailingOnly = TRUE)
file_gen <- here(str_remove(args[[1]], "^\\.\\./"))
file_smpl <- here(str_remove(args[[2]], "^\\.\\./"))

dir.create(path = here(args[[3]]), showWarnings = FALSE)
file_out_roh <- here(str_remove(args[[4]], "^\\.\\./"))
file_out_ftab <- here(str_remove(args[[5]], "^\\.\\./"))

gt_fmt <- args[[6]]     # [ "gp", "gt" ] note, that bcftools convert TAKES GT but EMITS GP format
rate_type <- args[[7]]  # [ "bp", "cM" ] check if your genotypes are living in "bp" or "cM" land
n_threads <- as.integer(args[[8]])

# --- running zooroh ---
# importing data
zoo_dat <- zoodata(genofile = file_gen,
                   supcol = 6,
                   poscol = 4,
                   chrcol = 1,
                   zformat = gt_fmt,
                   samplefile = file_smpl)

# for constant ROH segments class
default_rates <- c(2, 4, 8, 16, 32, 64, 128, 256, 512, 512) # use in CM land, last rate is non-HBD rate
bp_rates <- c(1,10,1e2,1e3,1e4,1e4) # use in bp land maybe?, last rate is non-HBD rate
rate_choices <- list(cM = default_rates, bp = bp_rates)

# set rate type
k_rates_used <- rate_choices[[rate_type]]
# define model
zoo_model <- zoomodel(K = length(k_rates_used), krates = k_rates_used)
# run
zoo_results <- zoorun(zoomodel = zoo_model, zooin = zoo_dat, nT = n_threads)

# --- exporting results ---

# zooplot_partitioning(list(test=zoo_results), ylim = c(0,1), nonhbd = FALSE)
# zooplot_hbdseg(zoo_results, chr = 1, coord=c(0e6,11.5e7))

# calculate inbreeding coefficients F_G−T
# relative to HBD classes and export as tsv
k_rates_used[1:(length(k_rates_used)-1)] |> 
  map(cumhbd, zres = zoo_results) |> 
  reduce(cbind) |> 
  cbind(zoo_dat@sample_ids) |> 
  as.data.frame() |> 
  set_names(nm = c(glue("F_{k_rates_used[1:(length(k_rates_used)-1)]}"), "sample_id")) |> 
  select(sample_id, everything()) |> 
  write_tsv(file_out_ftab)

# export identified HBD segments
rohbd(zoo_results) |> 
  as_tibble() |> 
  mutate(id = zoo_dat@sample_ids[id], 
         chrom = zoo_dat@chrnames[chrom]) |> 
  write_tsv(file_out_roh)
