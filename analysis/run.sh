#!/bin/bash

# Set the executable and command line parameters
PENNY=/path/to/GGSPenny 
GEOMETRY=/path/to/libGeoExample.so
DATACARD=TestRun.mac
OUTPUT=RunOutput.root

# Launch the simulation
time $PENNY -d $DATACARD -g $GEOMETRY -ro $OUTPUT