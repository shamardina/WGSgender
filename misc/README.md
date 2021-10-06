# Create counts for given reference genome

```sh
python create_nonN_counts.py /path/to/genome.fasta > genome.counts
cut -f 1,3 genome.counts > nonN_counts.csv
```

Only chromosomes 1-22, X and Y will be used for sex calling.
