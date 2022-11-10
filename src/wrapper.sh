#!/bin/bash

PROB=0.5
EDGES=4

for NODES in 10 12 15 18
do
  echo "$NODES ------------------------------------------------------"
  for SEED in $(seq 1 50)
  do
      echo "$SEED"
      python main.py "$NODES" "$EDGES" "$SEED" > log/watts_strogatz_"$NODES"_"$EDGES"_"$SEED".log
  done
done

