// --------------------------------------------------------------------------
//                            MC-LAMBDA INITIALIZE
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
//  Calls the initialize function LAMBDA_initialize, which is named for the 
//  alphabetically first entry-point function foo declared for code generation. 
//  Call the initialize function only once, even if you have multiple 
//  entry-point functions called in the function main.
//
// --------------------------------------------------------------------------
// Include Files
#include "mclambda_initialize.h"

// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//
// Arguments    : void
// Return Type  : void
//
// --------------------------------------------------------------------------
void mclambda_initialize()
{
  rt_InitInfAndNaN(8U);
}
