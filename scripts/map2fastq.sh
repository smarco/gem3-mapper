#!/bin/bash

awk '{print "@"$1"\n"$2"\n+\n"$3}' $1
