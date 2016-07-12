#!/bin/bash

scp Makefile* cluster2:/home/devel/smarco/workspace-C/gem-3.1.0
scp -r src/* cluster2:/home/devel/smarco/workspace-C/gem-3.1.0/src
scp -r include/* cluster2:/home/devel/smarco/workspace-C/gem-3.1.0/include
