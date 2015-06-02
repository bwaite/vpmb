#!/bin/bash

valgrind --leak-check=full build/vpmb vpm_decompression_input.json
