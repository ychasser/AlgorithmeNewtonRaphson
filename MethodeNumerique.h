#pragma once

#include "ModeleEANL.h"

class MethodesNumeriques
{

public : 

	static void newton_raph(int itemax, int* ite, double* crit_conver, double* crit_arret,
		double* sol, int derivee, double pas, int relax, enum STATUT* statut, ModeleEANL* modele, double* Fsol);


	static void MRINV(double **A, double **B, int N, int NRC, double &DETER, double EPS, double *X, int INDIC);


	static void JAC_NUM(double* x,double* f, double**A, ModeleEANL* modele,double pas);

	static double norme(double *x, int ndim);
	static double maxi(double* Tab, int ndim);

	static void ErreurMathLibC(int Numero);


};

