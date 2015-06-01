
# Building #

## C Version ##
   The C version requires that gcc, and make are installed, and CMake should
   also be installed. 

   From this folder you should run

   `cmake c/`

   to generate the Makefile. Then run

   `make`

   and everything should build properly.


## Python Version ##

   The python version doesn't need to be built, but does require python 2.6
   to be installed, since it needs the JSON library.

# Running #

## C ##
   The C program takes a json settings file as the first program argument. Optionally, 
   it can take a second argument, telling it to run the same way as the FORTRAN
   and Python program (this argument just has to exist, I was too lazy to check what it was).
   
   To run the simulation

   `./vpmb_c path_to_input_file`

   for the all at once

   `./vpmb_c path_to_input_file 1`

   The C program currently just outputs the step by step results of the simulation,
   and the final state of the dive_state structure. This is usually 5000 lines of output
   and up, so you probably want to pipe it to a file.

## Python ##

   The Python program has a few more options than the others. Most of the time, you just
   need to enter the directory, and run the program. It will load "vpm_decompression_input.json"
   and write the dive table to output.html. If you want to adjust the input or output files
   you can use the switches listed below.

   `-h, --help           show this help message and exit
   -i INPUT_FILE_NAME   Input file containing dive information
   -o OUTPUT_FILE_NAME  Output file for dive log
   -j JSON_OUTPUT       Output json instead of html`

## Testing ##

  You can run the current python test suite by going into the "tests" directory, and 
  running "runtests.py". If you're looking for sample input files to base a dive on,
  you can probably take something from one of the directories that's close to what
  you want, and modify it.

## Helper Scripts ##

  The "tools" directory contains a couple useful scripts.

  `json_to_fortran_file.py` takes a json settings file as input (-i flag) and generates 
  VPMDECO.IN, VPMDECO.SET, and ALTITUDE.SET files you can then provide to the 
  FORTRAN program.

  `validate_json.py` runs a number of checks on an input json file (again, -i flag),
  which can be useful for catching errors in your settings ahead of time.
