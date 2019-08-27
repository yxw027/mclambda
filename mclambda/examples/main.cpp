// --------------------------------------------------------------------------
//                         MAIN MC-LAMBDA SOFTWARE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//
//  DESCRIPTION:
//
//  This is the main software of use of the MC-LAMBDA routine. 
//  Its compilation allows to carry out a test as a test of the developed 
//  algorithm. (mclambda function). A default example is defined with predefined 
//  data a priori
//
//  It contains predefined functions for the creation of the input variables 
//  of the MC-LAMBDA routine, as well as a function that executes the entire 
//  creation and calculation process, which is called in the main function 
//  of the program
//
// --------------------------------------------------------------------------

// Include Files
#include <iostream>
#include "..\math_functions\rt_nonfinite.cpp"
#include "..\LAMBDA.cpp"
#include "main.h"
#include "..\LAMBDA_terminate.cpp"
#include "..\LAMBDA_emxAPI.cpp"
#include "..\LAMBDA_initialize.cpp"

using namespace std;

// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//
// Arguments    : double result[144]
// Return Type  : void
//
// --------------------------------------------------------------------------
static void argInit_12x12_real_T(double result[144])
{
  int idx0;
  int idx1;
  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < 12; idx0++) {
    for (idx1 = 0; idx1 < 12; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result[idx0 + 12 * idx1] = argInit_real_T();
    }
  }
}

// --------------------------------------------------------------------------
//
// Arguments    : double result[12]
// Return Type  : void
//
// --------------------------------------------------------------------------
static void argInit_12x1_real_T(double result[12])
{
  int idx0;
  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < 12; idx0++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx0] = argInit_real_T();
  }
}

// --------------------------------------------------------------------------
//
// Arguments    : void
// Return Type  : char
//
// --------------------------------------------------------------------------
static char argInit_char_T()
{
  return '?';
}

// --------------------------------------------------------------------------
//
// Arguments    : void
// Return Type  : double
//
// --------------------------------------------------------------------------
static double argInit_real_T()
{
  return 0.0;
}

// --------------------------------------------------------------------------
//
// Arguments    : void
// Return Type  : emxArray_char_T *
//
// --------------------------------------------------------------------------
static emxArray_char_T *c_argInit_UnboundedxUnbounded_c()
{
  emxArray_char_T *result;
  static int iv1[2] = { 2, 2 };

  int idx0;
  int idx1;
  // Set the size of the array.
  // Change this size to the value that the application requires.
  result = emxCreateND_char_T(2, *(int (*)[2])&iv1[0]);

  // Comment uncomment in order to change the method of calculation
  /*result->data[0] = 'P';
  result->data[1] = 'O';*/

  /*result->data[0] = 'M';
  result->data[1] = 'U';*/

  result->data[0] = 'N';
  result->data[1] = 'C';
  result->data[2] = 'A';
  result->data[3] = 'N';
  result->data[4] = 'D';
  result->data[5] = 'S';

  return result;
}

// --------------------------------------------------------------------------
//
// Arguments    : void
// Return Type  : void
//
// --------------------------------------------------------------------------
static void main_LAMBDA()
{
  // Definiton of variables and example values
  emxArray_real_T *afixed;
  emxArray_real_T *sqnorm;
  double ahat[12] = { -28490.8566886116, 65752.6299198198, 38830.3666554972, 5003.70833517778,
          -29196.0699104593, -297.658932458787, -22201.0284440701, 51235.8374755528,
          30257.7809603224, 3899.40332138829, -22749.1853575113, -159.278779870217 };
  double Qahat[144] = { 19068.8559508787, -15783.9722820370, -17334.2005875975, 14411.9239749603, 10055.7170089359, -14259.2952903872, 14858.8484050976, -12299.1993741839, -13507.1694819930, 11230.0704356810, 7835.62344938376, -11111.1393808147, 
          -15783.9722820370, 59027.7038409815, 38142.6927531102, 562.717388024645, -13830.0855960676, 27373.4263013019, -12299.1993747356, 45995.6129934030, 29721.5785731468, 438.480887460148, -10776.6902686912, 21329.9423774758,
          -17334.2005875975, 38142.6927531102, 28177.5653893528, -7000.50220497045, -11695.8674059306, 21886.1680630532, -13507.1694826246, 29721.5785738846, 21956.5440705992, -5454.93697674992, -9113.66310734779, 17054.1567378091,
          14411.9239749603, 562.717388024645, -7000.50220497045, 15605.5082283690, 5039.70281815470, -9648.96530646004, 11230.0704356773, 438.480887731461, -5454.93697653627, 12160.1358938811, 3927.04096307733, -7518.67445855756,
          10055.7170089359, -13830.0855960676, -11695.8674059306, 5039.70281815470, 6820.77250679480, -6880.24051213224, 7835.62344947055, -10776.6902682086, -9113.66310687634, 3927.04096320258, 5314.88728015545, -5361.22656658847,
          -14259.2952903872, 27373.4263013019, 21886.1680630532, -9648.96530646004, -6880.24051213224, 23246.5489626945, -11111.1393809211, 21329.9423779274, 17054.1567375591, -7518.67445829957, -5361.22656681708, 18114.1936088811,
          14858.8484050976, -12299.1993747356, -13507.1694826246, 11230.0704356773, 7835.62344947055, -11111.1393809211, 11578.3237340013, -9583.79156943782, -10525.0669778554, 8750.70438611838, 6105.68076067050, -8658.03053539344,
          -12299.1993741839, 45995.6129934030, 29721.5785738846, 438.480887731461, -10776.6902682086, 21329.9423779274, -9583.79156943782, 35840.7376978353, 23159.6717654859, 341.673569568934, -8397.42083743563, 16620.7344703582,
          -13507.1694819930, 29721.5785731468, 21956.5440705992, -5454.93697653627, -9113.66310687634, 17054.1567375591, -10525.0669778554, 23159.6717654859, 17108.9956804894, -4250.60009053988, -7101.55551676305, 13288.9534523001,
          11230.0704356810, 438.480887460148, -5454.93697674992, 12160.1358938811, 3927.04096320258, -7518.67445829957, 8750.70438611838, 341.673569568934, -4250.60009053988, 9475.43086798586, 3060.03207008500, -5858.70721928591,
          7835.62344938376, -10776.6902686912, -9113.66310734779, 3927.04096307733, 5314.88728015545, -5361.22656681708, 6105.68076067050, -8397.42083743563, -7101.55551676305, 3060.03207008500, 4141.47090961885, -4177.57899193454,
          -11111.1393808147, 21329.9423774758, 17054.1567378091, -7518.67445855756, -5361.22656658847, 18114.1936088811, -8658.03053539344, 16620.7344703582, 13288.9534523001, -5858.70721928591, -4177.57899193454, 14114.9563601479 };
  double method = 1.0;
  double param = 1;
  emxArray_char_T *type;
  double value = 10;
  double Ps;
  double Qzhat[144];
  double Z[144];
  double nfixed;
  double mu;
  emxInitArray_real_T(&afixed, 2);
  emxInitArray_real_T(&sqnorm, 2);

  // Initialize function input argument 'type'.
  type = c_argInit_UnboundedxUnbounded_c();

  // Call the entry-point 'mclambda'.
  mclambda(ahat, Qahat, method, param, type, value, afixed, sqnorm, &Ps, Qzhat, Z, &nfixed, &mu);

  // Print results of afixed vector
  for (int i=0; i<12; i++){
    cout<<afixed->data[i]<< endl;
  }
  emxDestroyArray_real_T(sqnorm);
  emxDestroyArray_real_T(afixed);
  emxDestroyArray_char_T(type);
}

// --------------------------------------------------------------------------
//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
// --------------------------------------------------------------------------
int main(int, const char * const [])
{
  // Initialize the application.
  // You do not need to do this more than one time.
  LAMBDA_initialize();

  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_LAMBDA();

  // Terminate the application.
  // You do not need to do this more than one time.
  LAMBDA_terminate();
  return 0;
}
// --------------------------------------------------------------------------
//                        End of mclambda software
// --------------------------------------------------------------------------