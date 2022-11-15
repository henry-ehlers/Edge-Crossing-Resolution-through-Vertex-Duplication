#!/bin/bash

# Iterate over first batch
for NODES in 8
do
  echo "$NODES --------------------------------------------------------------------------------------------------"
  for EDGES in 6
  do

    for SEED in $(seq 1 10)
    do
        echo "$SEED"
        python main.py "$NODES" "$EDGES" "$SEED" > log/watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
    done
  done
done