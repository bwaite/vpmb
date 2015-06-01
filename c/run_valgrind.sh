#!/bin/bash

valgrind --leak-check=full ../vpmb_c vpm_decompression_input.json
