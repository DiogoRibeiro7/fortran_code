#!/bin/bash

# This script compiles a Fortran source file using gfortran.
# Usage: ./compile_fortran.sh <input_file.f90> [accelerate]

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_file.f90> [accelerate]"
    exit 1
fi

input_file="$1"
output_name="${input_file%.f90}"
accelerate_option=""

if [ $# -gt 1 ] && [ "$2" == "accelerate" ]; then
    accelerate_option="-framework Accelerate"
fi

# Compile Fortran code using gfortran
gfortran -o "$output_name" "$input_file" $accelerate_option

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Executable '$output_name' created."
else
    echo "Compilation failed."
fi
