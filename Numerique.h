#pragma once
#include <iostream>
using namespace std;

double norme(double* x, int ndim);
double maxi(double* Tab, int ndim);
void JAC_NUM(double* x, double* f, double** jacob, int ndim, double pas, void (*syst) (double*, double*, double*), double* parametre = NULL);

void newton_raph(int dim, int itemax, int* ite, double* crit_conver, double* crit_arret,
	double* sol, int derivee, double pas, int relax,
	enum STATUT* statut, void (*syst) (double*, double*, double*),
	void (*JAC) (double*, double**, double*), double* Fsol, double* parametre = NULL);
