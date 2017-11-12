/* File : example.i */
%module GymSolver

%{
	#define SWIG_FILE_WITH_INIT
	#include <zlib.h>
	#include "GymSolver.h"
%}

// Get the NumPy typemaps
%include "numpy.i"

 // Get the STL typemaps
%include "stl.i"

// Handle standard exceptions
%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch (const std::out_of_range& e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
}
%init %{
  import_array();
%}

// Apply the 1D NumPy typemaps
%apply (float* INPLACE_ARRAY1, int DIM1) 
      {(float* array, int n)}
%apply (float* INPLACE_ARRAY1, int DIM1)
      {(float* array, int n), (float* pi, int m)}
%apply (float* INPLACE_ARRAY1, int DIM1)
      {(float* array, int n), (float* pi, int m), (float* v, int t)}
%apply (int DIM1  , float* INPLACE_ARRAY1)
      {(int length, float* data          )};
%apply (float** ARGOUTVIEW_ARRAY1, int* DIM1  )
      {(float** data             , int* length)};

/* Let's just grab the original header file here */
%include "GymSolver.h"


