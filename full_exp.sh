#!/bin/bash

set -e

TMP_DIR=tmp
RUN_CONFIDENCE_INTERVALS=0

mkdir -p tmp

pushd JTC
echo "Compiling the parameter estimation binary..."
./configure
make
popd

# Convert the CSV study data into the format that the parameter estimation
# binary expects.
echo -n "Preparing data from studies..."
python extract_data.py $TMP_DIR
echo "Done."

echo "Executing real data experiments."
rm -r tmp/real_data_runs

for i in {0..8}; do
  (python real_data_experiment.py $i && echo "Experiment $i done") &
done

wait $(jobs -p)

if [ "$RUN_CONFIDENCE_INTERVALS" = "1" ]
then
  echo "Running confidence intervals"

  for i in {0..8}; do
    python confidence_interval.py $i &
  done
fi
echo Done.
