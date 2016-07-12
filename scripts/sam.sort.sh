#!/bin/bash

awk '{if ($0 !~ /^@/) print}' $1 | sort -n
