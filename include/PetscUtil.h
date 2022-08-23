#ifndef PetscUtil_h
#define PetscUtil_h
#include <petsc.h>
#include "Home.h"

PetscErrorCode PetscInitialize(Home_t *home);
PetscErrorCode PetscTest(Home_t *home);

#endif