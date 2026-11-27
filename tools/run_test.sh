#!/bin/bash
# run_test.sh SOURCE INPUT NSHARDS  -> results/SOURCE/shard_<i>.tsv  (resumable; nice 19 is applied by the caller)
S=$1; IN=$2; K=$3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1  # one thread per process
cd "$(dirname "$0")/.."
mkdir -p results/$S
for i in $(seq 0 $((K-1))); do
  python3 -u tools/testall.py --n 8 --source $S --input $IN --out results/$S/shard_$i.tsv --shard $i --nshards $K &
done
wait
echo "DONE $S" >> results/$S/DONE
