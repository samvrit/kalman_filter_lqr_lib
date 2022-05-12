#!/bin/bash

echo "Compiling"

gcc -o obs observer_controller.c matrix_operations.c tests/test_observer_controller.c tests/project_specific.c -I. -Itests/ -DUNIT_TEST

echo "Running"

./obs.exe