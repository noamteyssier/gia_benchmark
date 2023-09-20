library(microbenchmark)
library(GenomicRanges)
library(dplyr)

load_bed <- function(filename) {
  x <- read.table(filename, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  with(x, GRanges(V1, IRanges(V2, V3)))
}

write_bed <- function(gr, filename) {
  write.table(as.data.frame(gr), file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}



closest <- function() {
  a.gr = load_bed(path_a)
  b.gr = load_bed(path_b)
  idx <- nearest(a.gr, b.gr)
  nearest <- b.gr[idx]
  write_bed(nearest, "/dev/null")
}

complement <- function() {
  a.gr = load_bed(path_large)
  ix <- gaps(a.gr)
  write_bed(ix, "/dev/null")
}

extend <- function() {
  a.gr = load_bed(path_large)
  ix <- resize(a.gr, width=1000)
  write_bed(ix, "/dev/null")
}

intersect <- function() {
  a.gr = load_bed(path_a)
  b.gr = load_bed(path_b)
  ix <- pintersect(a.gr, b.gr)
  write_bed(ix, "/dev/null")
}

merge <- function() {
  a.gr = load_bed(path_large)
  ix <- reduce(a.gr)
  write_bed(ix, "/dev/null")
}

run_sample <- function() {
  a.gr = load_bed(path_large)
  ix <- sample(a.gr, 100000)
  write_bed(ix, "/dev/null")
}

sort <- function() {
  a.gr = load_bed(path_large)
  sort.GenomicRanges(a.gr)
  write_bed(a.gr, "/dev/null")
}

subtract <- function() {
  a.gr = load_bed(path_a)
  b.gr = load_bed(path_b)
  ix <- setdiff(a.gr, b.gr)
  write_bed(ix, "/dev/null")
}

unit <- "ms"
times <- 30

path_a <- "~/projects/intervals/gia_bedtools_benchmark/data/a.bed"
path_b <- "~/projects/intervals/gia_bedtools_benchmark/data/b.bed"
path_large <- "~/projects/intervals/gia_bedtools_benchmark/data/large.bed"

bench_closest <- microbenchmark(closest(), times=times, unit=unit)
bench_complement <- microbenchmark(complement(), times=times, unit=unit)
bench_extend <- microbenchmark(extend(), times=times, unit=unit)
bench_intersect <- microbenchmark(intersect(), times=times, unit=unit)
bench_merge <- microbenchmark(merge(), times=times, unit=unit)
bench_sample <- microbenchmark(run_sample(), times=times, unit=unit)
bench_sort <- microbenchmark(sort(), times=times, unit=unit)
bench_subtract <- microbenchmark(subtract(), times=times, unit=unit)

timing <- bind_rows(
  data.frame(bench_closest),
  data_frame(bench_complement),
  data_frame(bench_extend),
  data_frame(bench_intersect),
  data_frame(bench_merge),
  data_frame(bench_sample),
  data_frame(bench_sort),
  data_frame(bench_subtract),
) %>% 
  mutate(
    expr = gsub("\\(\\)", "", expr),
    time = time / 1e9
  ) %>% 
  mutate(
    expr = gsub("run_", "", expr)
  )

timing

agg <- timing %>% 
  group_by(expr) %>% 
  summarize(
    mean = mean(time),
    stddev = sd(time),
    median = median(time),
    min = min(time),
    max = max(time),
  ) %>% 
  rename(
    category = expr
  ) %>% 
  mutate (
    tool = "granges",
    named = "named",
    handle = "granges",
    user = NA,
    system = NA,
  )

agg

write.csv(
  agg, 
  file = "~/projects/intervals/gia_bedtools_benchmark/results/granges.csv",
  row.names = F
  )
