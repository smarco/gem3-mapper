#!/bin/bash

awk '{if ($5=="-") print}' $1
