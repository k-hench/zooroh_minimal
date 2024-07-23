convert to actual rescale map (wrong_ref):

```sh
zgrep mscaf_a1_01 normalized_genetic_map.tsv.gz | \
  cut -f 1-4,10 | \
  sed 's/mscaf_a1_01/NC_072356.1/g' | \
  gzip > genetic_map_wrongref_norm.bed.gz 
```

translate bp to cM

```sh
py/bp2cm \
  -i ../results/genotypes/mirang_test.gen.gz \
  -m ../data/genetic_map_wrongref_norm.bed.gz \
  -t 4 \
  -o ../results/genotypes/mirang_test_cM.gen.gz
```