#!/bin/bash

scp -P 10000 Makefile* aopccuda.uab.es:/home/smarco/gem3-gpu.Dev
scp -P 10000 -r src/* aopccuda.uab.es:/home/smarco/gem3-gpu.Dev/src
scp -P 10000 -r include/* aopccuda.uab.es:/home/smarco/gem3-gpu.Dev/include
