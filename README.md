# gia_benchmark

a benchmark of gia and other gene interval tools on various functionality.

## Information

This was originally done using the following versions:

```text
gia == 0.1.16
bedtools == 2.31.0
bedops == 2.4.41
GenomicRanges == 3.17
```

## Running the Benchmark

### Get-Fasta

To recreate the get-fasta command you'll need to download the homo-sapiens primary
assembly:

You can download it easily with [`ggetrs`](https://noamteyssier.github.io/ggetrs/ensembl/ref.html)

``` bash
cd data/
ggetrs ref -Dd dna
```

and generate its indexed fasta (.fai)

```bash
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cd ../
```

### Benchmark

The benchmarks are run using [`just`](https://just.systems/) and [`hyperfine`](https://github.com/sharkdp/hyperfine).

You'll need to install each of those:

```bash
cargo install just hyperfine
```

Then you can just run:

```bash
just bench
```

And all the benchmarks should run.

## Analysis

The visualization of the data can be viewed using the provided notebooks
which recreate the analysis within the manuscript.
