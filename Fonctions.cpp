/*--------------------------------------------------------------------------*/
/* Définition des fonctions définissant les systèmes à résoudre             */
/*   et leurs jacobiens analytiques				                            */
/*--------------------------------------------------------------------------*/

#include <math.h>
#include <iostream>
using namespace std;

// Système de test 1
void syst1( double *x, double *F1, double *parametre)
{
	/* Définition de la fonction */
	F1[1]=x[1]*x[1]-x[1]*x[2]*x[2]-2;
	F1[2]=2*x[1]*x[1]-3*x[1]*x[2]*x[2]+3;
}

// Jacobien de test 1
void JAC1(double* x, double** jacob, double* parametre)
{
	/* Définition du jacobien */
	jacob[1][1] = 2 * x[1] - pow(x[2], 2);
	jacob[1][2] = -2 * x[1];
	jacob[2][1] = 4 * x[1] - 3 * pow(x[2], 2);
	jacob[2][2] = -6 * x[1] * x[2];
}

// Système de test 2
void syst2(double* x, double* F1, double* parametre)
{
	/* Définition de la fonction */
	F1[1] = cos(x[1]) - x[1] * x[1] * x[1];
}

// Jacobien de test 2
void JAC2(double* x, double** jacob, double* parametre)
{
	/* Définition du jacobien */
	jacob[1][1] = -sin(x[1]) - 3 * x[1] * x[1];
}



