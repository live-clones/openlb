#!/usr/bin/env bash

COUNTER=0
COUNTER_FAILED=0
COUNTER_SUCCESS=0

for target_file in expressions/*.out; do
    stripped_name="${target_file#expressions/}" # Remove "expressions/" prefix
    base_name="${stripped_name%.out}" # Remove ".txt" suffix
    out_file="expressions/include/$base_name.cse.h" # Rebuild the path for obj_file
    let COUNTER++
    if [ ! -f "$out_file" ]; then
        echo "Failed to optimize: $target_file"
        let COUNTER_FAILED++
    fi
done

declare -i COUNTER_SUCCESS=COUNTER-COUNTER_FAILED 
echo "Summary: SUCCESS=$COUNTER_SUCCESS / FAILED=$COUNTER_FAILED / TOTAL=$COUNTER"
