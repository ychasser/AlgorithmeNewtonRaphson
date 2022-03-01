#include <iostream>
#include "MathLibC.h"
#include "Numerique.h"
#include "type.h"
#include "math.h"

using namespace std;

/*--------------------------------------------------------------------------*/
/* Définition du jacobien numérique										*/
/*--------------------------------------------------------------------------*/
void JAC_NUM(double* x, double* f, double** jacob, int ndim, double pas,
	void (*syst) (double*, double*, double*), double* parametre)
{
	int i, j;
	double* g;
	double u;

	// Allocation dynamique de g
	// -------------------------
	g = new double[ndim + 1];

	for (j = 1; j <= ndim; j++)
	{
		u = x[j];
		x[j] = x[j] + x[j] * pas;

		for (i = 1; i <= ndim; i++)
		{
			syst(x, g, parametre);
			jacob[i][j] = (g[i] - f[i]) / pas / u;
		}
		x[j] = u;
	}
	delete[] g;
}


/*-------------------------------------------------------------------------------*/
/*			Fonction permettant de calculer la norme Euclidienne				 */
/*					d'un vecteur de dimension ndim								 */
/*-------------------------------------------------------------------------------*/
double norme(double *x, int ndim)
{
	int i;
	double norm = 0;

	for (i=1; i<=ndim; i++)
	{
		norm = norm +  pow(x[i],2);
	}

	norm = sqrt (norm);
	return (norm);
}

/*----------------------------------------------------------------------------*/
/* Fonction visant à rechercher la composante maximal d'un vecteur de réels  */
/*----------------------------------------------------------------------------*/
double maxi(double* Tab, int ndim)
{
	int i;
	double max;

	max = Tab[1];

	for (i = 2; i <= ndim; i++)
	{
		if (Tab[i] > max)
		{
			max = Tab[i];
		}
	}
	return (max);
}


/*--------------------------------------------------------------------------*/
/* Sous-programme de résolution d'un système d'équations non-linéaires		*/
/* par la méthode de Newton-Raphson											*/
/*																			*/
/* AUTEURS :G.HETREUX 								               			*/
/* DATE : 25/02/2001														*/
// ARGUMENTS																//
//																			//
// int dim : Dimension du système à résoudre                                //
// int itemax : Nombre maximum d'itérations si                              //
//              aucune convergence/stagnation                               //
//                                                                          //
// int* ite : compteur d'iterations                                         //
// double* crit_conver : critère de convergence (epsilon 1)                 //
// double* crit_arret : critère de stagnation (epsilon 2)                   //
// double* sol : vecteur solution calculé au cours de la méthode            //
// int derivee : choix sur le mode de calcul                                //
//               du jacobien (1 : Analytique, 2: Numérique)                 //
//                                                                          //
// double pas : pas utilisé pour la calcul du jacobien numérique            //
// int relax : choix sur le mode de relaxation                              //
//             (1 : pas de relaxation, 2 : relaxation numérique)            // 
//                                                                          //
// enum STATUT* statut : Indicateur sur le statut de la méthode             //
//                       (EN_COURS, CONVERGENCE,                            //
//                       STAGNATION_DU_PAS, NBRE_ITERATION_DEPASSE)         //
//                                                                          //
// void (*syst) (double*, double*, double*) : Nom de la fonction            //
//                                            représentant le système à     //
//                                            résoudre (voir Fonctions.cpp) //
//                                                                          //
// void (*syst) (double*, double*, double*) : Nom de la fonction            //
//                         représentant le jacobien analytique du système   //
//                         à résoudre (voir Fonctions.cpp)                  //
//                                                                          //
// double* Fsol : Valeur du système à la solution                           //
// double* parametre : vecteur de paramètres à passer aux systèmes          //
//                     paramétrés                                           //
//                                                                          //
// VARIABLES                                                                //
//                                                                          //
// int i : compteur                                                         //
// Deter : Indique si le déterminant de la matrice                          //
//         jacobienne a pu etre calculé                                     //
//                                                                          //
// Norme : Norme de Fsol (ou F_relax)                                       //
// alpha : coefficient de la relaxation numérique                           //
// pivot_min : zero de la méthode du pivot (DMRINV)                         //
// type_calcul : Précision sur le choix de calcul                           //
//               à employer pour MRINV (voir MathLibC.h)                    //
//                                                                          //
// A : Matrice à inverser représentant le système d'equations linéaires     //
// H : Pas calculé par résolution du système d'équations linéaires          //
// Y : Vecteur solution mis à jour                                          //
// Evol : Vecteur représentant l'évolution de la solution                   //
//        d'une itération sur l'autre                                       //
//                                                                          //
// F_relax : Valeur du système au cours de la relaxation                    //
//                                                                          //
/*--------------------------------------------------------------------------*/
void newton_raph(int ndim, int itemax, int* ite, double* crit_conver, double* crit_arret,
	double* sol, int derivee, double pas, int relax, enum STATUT* statut, void (*syst) (double*, double*, double*),
	void (*JAC) (double*, double**, double*), double* Fsol, double* parametre)
{
	int i;
	double Deter;
	double Norme;
	double alpha;

	const double pivot_min = 1e-30;
	const int type_calcul = 1;

	// Declaration et allocation dynamique des tableaux de travail 
	double** A; // matrice de dimension (n, n+1)
	double* H; // vecteur de dimension n
	double* Y; // vecteur de dimension n
	double* Evol; // vecteur de dimension 
	double* F_relax;

	H = new double[ndim + 1];
	Y = new double[ndim + 1];
	Evol = new double[ndim + 1];
	F_relax = new double[ndim + 1];
	
	A = new double* [ndim + 1];
	for (i = 1; i <= ndim; i++)
	{
		A[i] = new double[ndim + 1 + 1];
	}


	// Début de la boucle de calcul - Initialisation
	*ite = 0;
	syst(sol, Fsol, parametre); // Intitialisation de F
	Norme = norme(Fsol, ndim); // Initialisation de la norme de F
	*statut = EN_COURS; // Initialisation du statut de la méthode

	do
	{
		//  Test sur le critère de convergence
		if (Norme > *crit_conver)
		{
			(*ite)++; // Incrément du nombre d'itérations
			if (*ite <= itemax)
			{
				// Calcul du Jacobien
				switch (derivee)
				{
				case 1: // Jacobien analytique
					JAC(sol, A, parametre);
					break;
				case 2: // Jacobien numérique
					JAC_NUM(sol, Fsol, A, ndim, pas, syst, parametre);
					break;
				}
				// Remplissage de la matrice A : définition du système d'équations linéaires
				for (i = 1; i <= ndim; i++)
				{
					A[i][ndim + 1] = -Fsol[i];
				}

				// Appel a Mrinv : procédure de résolution de systeme d'equations lineaires
				MRINV(A, A, ndim, 0, Deter, pivot_min, H, type_calcul);

				// Calcul de sol au pas suivant
				switch (relax)
				{
				case 1: // Pas de ralaxation
					for (i = 1; i <= ndim; i++)
					{
						Y[i] = sol[i] + H[i];
					}
					break;
				case 2: // Relaxation numérique
					alpha = 1;
					do {
						for (i = 1; i <= ndim; i++)
						{
							Y[i] = sol[i] + alpha * H[i];
						}
						syst(Y, F_relax, parametre);
						alpha = alpha / 2.;
					} while (norme(F_relax, ndim) > Norme);

				}


				// Calcul de l'évolution entre deux itérations
				for (i = 1; i <= ndim; i++)
				{
					Evol[i] = fabs((Y[i] - sol[i]) / sol[i]);
				}

				// Si la méthode ne stagne pas, on met à jour la solution
				if (maxi(Evol, ndim) > *crit_arret)
				{
					for (i = 1; i <= ndim; i++)
					{
						sol[i] = Y[i];
					}
					syst(sol, Fsol, parametre); // Mise à jour de F
					Norme = norme(Fsol, ndim); // Mise à jour de la norme de F
					cout << Norme << endl;
				}
				else
				{
					*statut = STAGNATION_DU_PAS;
				}
			}
			else
			{
				*statut = NBRE_ITERATION_DEPASSE;
			}

		}
		else
		{
			*statut = CONVERGENCE;
		}

	} 	while (*statut == EN_COURS);

	// Mise à jour des paramètres de sortie
	*crit_arret = maxi(Evol, ndim);
	*crit_conver = Norme;

	// Destruction des structures dynamiques
	delete[] H;
	delete[] Evol;
	delete[] Y;
	delete[] F_relax;

	for (i = 1; i <= ndim; i++)
	{
		delete[] A[i];
	}
	delete[] A;
}
