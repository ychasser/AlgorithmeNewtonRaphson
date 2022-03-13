#pragma once
#include "ModeleEANL.h"

class ModeleSysteme1 : public ModeleEANL
{
public : 

	ModeleSysteme1(int ndim);

	void EvaluerF(double* X,double* F) ;

	void EvaluerJacob(double* X,double** Jacob) ;



};

