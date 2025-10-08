#!/bin/sh
./cavity3d  --RESOLUTION 128 --TIME_STEPS 100 --NO_EXPORT_RESULTS 1 2>&1 | tee output.log
