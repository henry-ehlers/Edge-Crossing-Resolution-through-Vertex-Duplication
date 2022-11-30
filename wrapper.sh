#!/bin/bash

# Iterate over first batch
for NODES in 20
do
  echo "$NODES --------------------------------------------------------------------------------------------------"
  for EDGES in 8
  do

    for SEED in $(seq 1 10)
    do
        echo "$SEED"
        python main_inner.py "$NODES" "$EDGES" "$SEED" > log/inner_watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
        python main_outer.py "$NODES" "$EDGES" "$SEED" > log/outer_watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
    done
  done
done