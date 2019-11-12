#!/bin/bash

rm -rf results_0 results_1
mpirun \
-np 10 ./hemepure0 -in input_0.xml -out results_0 : \
-np 10 ./hemepure1 -in input_1.xml -out results_1
