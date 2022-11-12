#!/bin/bash

# Iterate over first batch
for NODES in 12 13
do
  echo "$NODES --------------------------------------------------------------------------------------------------"
  for EDGES in 6 8
  do

    for SEED in $(seq 51 100)
    do
        echo "$SEED"
        python main.py "$NODES" "$EDGES" "$SEED" > log/watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
    done
  done
done