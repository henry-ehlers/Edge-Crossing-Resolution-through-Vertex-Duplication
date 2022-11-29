#!/bin/bash

# Iterate over first batch
for NODES in 15
do
  echo "$NODES --------------------------------------------------------------------------------------------------"
  for EDGES in 4, 6
  do

    for SEED in $(seq 1 20)
    do
        echo "$SEED"
        python main_inner.py "$NODES" "$EDGES" "$SEED" > log/watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
    done
  done
done