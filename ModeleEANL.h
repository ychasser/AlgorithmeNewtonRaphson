#pragma once
#ifndef MODELE_EANL_H
#define MODELE_EANL_H

class ModeleEANL
{

protected:

	int dimension;
	int dimParamEntier;
	int dimParamReel;

	double* parametreInt;

	double* parametreReel;


public:

	ModeleEANL(int dim);

	virtual void EvaluerF(double* X,double* F) = 0;

	virtual void EvaluerJacob(double* X, double** Jacob) = 0;

	int Get_dimension();



};

#endif




