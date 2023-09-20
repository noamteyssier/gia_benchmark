
BED_A := "data/a.bed"
BED_B := "data/b.bed"
SORTED_A := "data/a.sorted.bed"
SORTED_B := "data/b.sorted.bed"
BED_LARGE := "data/large.bed"
GENOME := "data/genome.txt"
FASTA := "data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

MAX_CHROMOSOME_LEN := "10000000"
NUM_INTERVALS := "200000"
NUM_INTERVALS_LARGE := "500000"
NUM_CHROMOSOMES := "10"
INTERVAL_SIZE := "1000"
EXTEND_SIZE := "1000"
NUM_SAMPLES := "100000"

RESULTS_DIR := "results"



bench: bench_intersect bench_merge bench_sample bench_extend bench_random bench_subtract bench_sort bench_closest bench_complement bench_getfasta


gen_data:
  gia random -n {{NUM_INTERVALS}} -l 250 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -o {{BED_A}} -s 0
  gia random -n {{NUM_INTERVALS}} -l 500 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -o {{BED_B}} -s 1

gen_sorted_data:
  gia random -n {{NUM_INTERVALS}} -l 250 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -s 3 | \
    gia sort -o {{SORTED_A}}
  gia random -n {{NUM_INTERVALS}} -l 500 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -s 4 | \
    gia sort -o {{SORTED_B}}

gen_data_large:
  gia random -n {{NUM_INTERVALS_LARGE}} -l 250 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -s 2 | \
  gia sort -o {{BED_LARGE}}

gen_genome: gen_data_large
  for i in $(seq 1 {{NUM_CHROMOSOMES}}); do \
    echo $i $(echo "{{MAX_CHROMOSOME_LEN}}+{{EXTEND_SIZE}}+1" | bc); \
  done | \
    tr " " "\t" > {{GENOME}}  

gen_sorted_data_matrix:
  for i in 100 250 500 1000 2500 5000 10000 25000 50000 100000 250000 500000 750000 1000000 2500000 5000000 10000000; do \
    gia random -n $i -l 250 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -s 0 | \
      gia sort -o data/random_a_${i}.bed; \
    gia random -n $i -l 500 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -s 1 | \
      gia sort -o data/random_b_${i}.bed; \
  done
  
gen_sorted_data_matrix_high:
  for i in 10 25 50 75; do \
    for j in 10 100 1000 10000 100000 1000000; do \
      num=$(echo $i*$j | bc); \
      gia random -n $num -l 250 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -s 0 | \
        gia sort -o data/random_a_${num}.bed; \
      gia random -n $num -l 500 -c {{NUM_CHROMOSOMES}} -m {{MAX_CHROMOSOME_LEN}} -s 1 | \
        gia sort -o data/random_b_${num}.bed; \
    done \
  done

#############
# INTERSECT #
#############

bench_intersect: gen_data bench_intersect_fractional bench_intersect_sorted

bench_intersect_sorted: gen_sorted_data
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/intersect.csv" \
    "gia intersect -a {{SORTED_A}} -b {{SORTED_B}} -o /dev/null" \
    "gia intersect -N -a {{SORTED_A}} -b {{SORTED_B}} -o /dev/null" \
    "bedtools intersect -a {{SORTED_A}} -b {{SORTED_B}} > /dev/null" \
    "bedops -i {{SORTED_A}} {{SORTED_B}} > /dev/null";
  

bench_intersect_fractional: gen_data
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/insersect_fractional.csv" \
    "gia intersect -f 0.8 -a {{BED_A}} -b {{BED_B}} -o /dev/null" \
    "gia intersect -f 0.8 -N -a {{BED_A}} -b {{BED_B}} -o /dev/null" \
    "bedtools intersect -f 0.8 -a {{BED_A}} -b {{BED_B}} > /dev/null";

bench_intersect_query_unique: gen_data
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/insersect_query_unique.csv" \
    "gia intersect -qu -a {{BED_A}} -b {{BED_B}} -o /dev/null" \
    "gia intersect -qu -N -a {{BED_A}} -b {{BED_B}} -o /dev/null" \
    "bedtools intersect -wa -u -a {{BED_A}} -b {{BED_B}} > /dev/null";

############
# SUBTRACT #
############

bench_subtract: gen_sorted_data
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/subtract.csv" \
    "gia subtract -a {{SORTED_A}} -b {{SORTED_B}} -o /dev/null" \
    "gia subtract -N -a {{SORTED_A}} -b {{SORTED_B}} -o /dev/null" \
    "bedtools subtract -a {{SORTED_A}} -b {{SORTED_B}} > /dev/null" \
    "bedops -d {{SORTED_A}} {{SORTED_B}} > /dev/null";

#############
#   MERGE   #
#############

bench_merge: gen_data_large
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/merge.csv" \
    "gia merge -i {{BED_LARGE}} -o /dev/null" \
    "gia merge -N -i {{BED_LARGE}} -o /dev/null" \
    "bedtools merge -i {{BED_LARGE}} > /dev/null" \
    "bedops -m {{BED_LARGE}} > /dev/null";


############
# SAMPLING #
############

bench_sample: gen_data_large
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/sample.csv" \
    "gia sample -i {{BED_LARGE}} -n {{NUM_SAMPLES}} -o /dev/null" \
    "gia sample -N -i {{BED_LARGE}} -n {{NUM_SAMPLES}} -o /dev/null" \
    "bedtools sample -i {{BED_LARGE}} -n {{NUM_SAMPLES}} > /dev/null";

##########
# EXTEND #
##########

bench_extend: gen_data_large gen_genome
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/extend.csv" \
    "gia extend -i {{BED_LARGE}} -b {{EXTEND_SIZE}} -o /dev/null" \
    "gia extend -N -i {{BED_LARGE}} -b {{EXTEND_SIZE}} -o /dev/null" \
    "bedtools slop -i {{BED_LARGE}} -g {{GENOME}} -b {{EXTEND_SIZE}} > /dev/null" \
    "bedops --range {{EXTEND_SIZE}} -u {{BED_LARGE}} > /dev/null";


##########
# RANDOM #
##########

bench_random: gen_genome
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/random.csv" \
    "gia random -n {{NUM_INTERVALS_LARGE}} -l {{INTERVAL_SIZE}} -g {{GENOME}} -o /dev/null" \
    "bedtools random -l {{INTERVAL_SIZE}} -n {{NUM_INTERVALS_LARGE}} -g {{GENOME}} > /dev/null";

########
# SORT #
########

bench_sort: gen_data_large
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/sort.csv" \
    "gia sort -i {{BED_LARGE}} -o /dev/null" \
    "gia sort -N -i {{BED_LARGE}} -o /dev/null" \
    "bedtools sort -i {{BED_LARGE}} > /dev/null" \
    "sort-bed {{BED_LARGE}} > /dev/null";

###########
# CLOSEST #
###########

bench_closest: gen_sorted_data
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/closest.csv" \
    "gia closest -a {{SORTED_A}} -b {{SORTED_B}} -o /dev/null" \
    "gia closest -N -a {{SORTED_A}} -b {{SORTED_B}} -o /dev/null" \
    "bedtools closest -a {{SORTED_A}} -b {{SORTED_B}} -t first > /dev/null" \
    "closest-features --closest {{SORTED_A}} {{SORTED_B}} > /dev/null";


##############
# COMPLEMENT #
##############

bench_complement: gen_data_large gen_genome
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/complement.csv" \
    "gia complement -i {{BED_LARGE}} -o /dev/null" \
    "gia complement -N -i {{BED_LARGE}} -o /dev/null" \
    "bedtools complement -i {{BED_LARGE}} -g {{GENOME}} > /dev/null" \
    "bedops -c {{BED_LARGE}} > /dev/null";

#############
# GET FASTA #
#############

bench_getfasta: gen_sorted_data
  hyperfine \
    -u millisecond \
    --warmup 10 \
    --export-csv "{{RESULTS_DIR}}/getfasta.csv" \
    "gia get-fasta -b {{SORTED_A}} -f {{FASTA}} -o /dev/null" \
    "bedtools getfasta -fi {{FASTA}} -bed {{SORTED_A}} > /dev/null";

###################
# INTERSECT RANGE #
###################

bench_intersect_range: gen_sorted_data_matrix
  for i in 100 250 500 1000 2500 5000 10000 25000 50000 100000 250000 500000 750000 1000000 2500000 5000000 10000000; do \
    hyperfine \
      --warmup 1 \
      -r 3 \
      -u millisecond \
      --export-csv "{{RESULTS_DIR}}/intersect_range_${i}.csv" \
      "gia intersect -a "data/random_a_${i}.bed" -b "data/random_b_${i}.bed" -o /dev/null" \
      "gia intersect -S -a "data/random_a_${i}.bed" -b "data/random_b_${i}.bed" -o /dev/null" \
      "bedtools intersect -a "data/random_a_${i}.bed" -b "data/random_b_${i}.bed" > /dev/null" \
      "bedtools intersect -a "data/random_a_${i}.bed" -b "data/random_b_${i}.bed" -sorted > /dev/null" \
      "bedops -i "data/random_a_${i}.bed" "data/random_b_${i}.bed" > /dev/null"; \
  done

bench_intersect_high_range: gen_sorted_data_matrix_high
  for i in 10 25 50 75; do \
    for j in 10 100 1000 10000 100000 1000000; do \
      num=$(echo $i*$j | bc); \
      hyperfine \
        --warmup 1 \
        -r 3 \
        -u millisecond \
        --export-csv "{{RESULTS_DIR}}/intersect_range_high_${num}.csv" \
        "gia intersect -S -a "data/random_a_${num}.bed" -b "data/random_b_${num}.bed" -o /dev/null" \
        "bedops -i "data/random_a_${num}.bed" "data/random_b_${num}.bed" > /dev/null"; \
    done \
  done
