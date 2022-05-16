#pragma once

#include "ModeleEANL.h"
#include "ModeleSysteme1.h"

class MethodesNumeriques
{

public : 

	// Declaration des methodes NawtonRaphson, JacobienNumerique et Norme

	static void MRINV(double **A, double **B, int N, int NRC, double &DETER, double EPS, double *X, int INDIC);

	static double maxi(double* Tab, int ndim);

	static void ErreurMathLibC(int Numero);

	static void newton_raph(int itemax, int* ite, double* crit_conver, double* crit_arret, double* X, int derivee, double pas, int relax, STATUT* statut, ModeleEANL* modele, double* F);

	static double norme(double* x, int ndim);

	static void JAC_NUM(double* x, double* f, double** jacob, ModeleEANL* modele, double pas);

};

