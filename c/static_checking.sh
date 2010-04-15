#!/bin/bash

cppcheck vpmb.c
rats vpmb.c
splint vpmb.c
~/Jars/pmd-4.2.5/bin/cpd.sh .