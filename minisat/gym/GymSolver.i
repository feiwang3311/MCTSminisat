/* File : example.i */
%module GymSolver

%{
#include <zlib.h>
#include "GymSolver.h"
%}

/* Let's just grab the original header file here */
%include "GymSolver.h"
