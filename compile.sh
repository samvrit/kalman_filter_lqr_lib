#!/bin/bash

echo "Cleaning"

rm obs.exe

echo "Compiling"

gcc -o obs observer_controller.c matrix_operations.c tests/test_observer_controller.c -I. -Itests/ -DUNIT_TEST

echo "Running"

./obs.exe