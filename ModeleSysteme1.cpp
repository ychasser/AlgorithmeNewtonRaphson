#include <cmath>
#include "ModeleSysteme1.h"


ModeleSysteme1::ModeleSysteme1(int ndim) : ModeleEANL(ndim)
{
	
}



void ModeleSysteme1::EvaluerF(double* x,double* valueF)
{
	/* D�finition de la fonction */
	// @TODO : coder ici les equations
	valueF[1] = x[1] * x[1] - x[1] * x[2] * x[2] - 2;
	valueF[2] = 2 * x[1] * x[1] - 3 * x[1] * x[2] * x[2] + 3;
}

void ModeleSysteme1::EvaluerJacob(double* x,double** valueJacob)
{
	/* D�finition du jacobien */
	// @TODO : coder ici le jacobien
	valueJacob[1][1] = 2 * x[1] - pow(x[2], 2);
	valueJacob[1][2] = -2 * x[1];
	valueJacob[2][1] = 4 * x[1] - 3 * pow(x[2], 2);
	valueJacob[2][2] = -6 * x[1] * x[2];
}

