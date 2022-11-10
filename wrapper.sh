#!/bin/bash

# Specify Hyper-arameters
PROB=0.5

# Specify new edge attachment parameter
EDGES=4

# Iterate over first batch
for NODES in 5 8 10 12 15
do
  echo "$NODES --------------------------------------------------------------------------------------------------"
  for SEED in $(seq 1 50)
  do
      echo "$SEED"
      python main.py "$NODES" "$EDGES" "$SEED" > log/watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
  done
done

# Specify new edge attachment parameter
EDGES=6

# Iterate over first batch
for NODES in 8 10 12 15
do
  echo "$NODES --------------------------------------------------------------------------------------------------"
  for SEED in $(seq 1 50)
  do
      echo "$SEED"
      python main.py "$NODES" "$EDGES" "$SEED" > log/watts_strogatz_"$NODES"_"$EDGES"_"$SEED"_"$PROB".log
  done
done
