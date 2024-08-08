#!/bin/bash

# Run clang-format on all C++ source and header files except those in
# the build directory

if [ $# -eq 0 ]; then
    echo "Usage: bash clang_format_check.sh <clang-format-executable>"
    exit 1
fi

CLANG=$1

if [ ! -f $(which $CLANG) ]; then
    echo "$CLANG could not be found"
    exit 1
fi

output=$(find -type d -name build -prune -false -o \
              -regextype posix-extended -regex ".*\.(cpp|h)" \
              -exec sh -c '$1 -style=file $0 | diff -u $0 -' {} $CLANG \;)

if [[ $output ]]; then
    echo "$output"
    echo "Error: source code formatting errors ($($CLANG --version))"
    exit 1
else
    echo "Success: no source code formatting errors ($($CLANG --version))"
fi

exit 0
