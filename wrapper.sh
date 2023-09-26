#!/bin/bash

# Iterate over first batch
for NODES in 14
do
  echo "$NODES --------------------------------------------------------------------------------------------------"
  for EDGES in 1
  do

    for SEED in $(seq 90 100)
    do
        echo "$NODES $SEED"
        python main_inner.py "$NODES" "$EDGES" "$SEED" > log/inner_complete_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
        # python main_outer.py "$NODES" "$EDGES" "$SEED" > log/outer_watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
    done
  done
done