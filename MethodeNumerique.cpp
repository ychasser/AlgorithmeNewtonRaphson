#include <cmath>
#include <iostream>
using namespace std;
#include "type.h"
#include "MethodeNumerique.h"


/*--------------------------------------------------------------------------*/
/* Sous-programme de r�solution d'un syst�me d'�quations non-lin�aires		*/
/* par la m�thode de Newton-Raphson											*/
/*																			*/
/* AUTEURS :G.HETREUX 								               			*/
/* DATE : 25/02/2001														*/
// ARGUMENTS																//
//																			//
// int dim : Dimension du syst�me � r�soudre                                //
// int itemax : Nombre maximum d'it�rations si                              //
//              aucune convergence/stagnation                               //
//                                                                          //
// int* ite : compteur d'iterations                                         //
// double* crit_conver : crit�re de convergence (epsilon 1)                 //
// double* crit_arret : crit�re de stagnation (epsilon 2)                   //
// double* sol : vecteur solution calcul� au cours de la m�thode            //
// int derivee : choix sur le mode de calcul                                //
//               du jacobien (1 : Analytique, 2: Num�rique)                 //
//                                                                          //
// double pas : pas utilis� pour la calcul du jacobien num�rique            //
// int relax : choix sur le mode de relaxation                              //
//             (1 : pas de relaxation, 2 : relaxation num�rique)            // 
//                                                                          //
// enum STATUT* statut : Indicateur sur le statut de la m�thode             //
//                       (EN_COURS, CONVERGENCE,                            //
//                       STAGNATION_DU_PAS, NBRE_ITERATION_DEPASSE)         //
//                                                                          //
// void (*syst) (double*, double*, double*) : Nom de la fonction            //
//                                            repr�sentant le syst�me �     //
//                                            r�soudre (voir Fonctions.cpp) //
//                                                                          //
// void (*syst) (double*, double*, double*) : Nom de la fonction            //
//                         repr�sentant le jacobien analytique du syst�me   //
//                         � r�soudre (voir Fonctions.cpp)                  //
//                                                                          //
// double* Fsol : Valeur du syst�me � la solution                           //
// double* parametre : vecteur de param�tres � passer aux syst�mes          //
//                     param�tr�s                                           //
//                                                                          //
// VARIABLES                                                                //
//                                                                          //
// int i, ite_relax : compteurs                                             //
// Deter : Indique si le d�terminant de la matrice                          //
//         jacobienne a pu etre calcul�                                     //
//                                                                          //
// Norme : Norme de Fsol (ou F_relax)                                       //
// alpha : coefficient de la relaxation num�rique                           //
// pivot_min : zero de la m�thode du pivot (DMRINV)                         //
// type_calcul : Pr�cision sur le choix de calcul                           //
//               � employer pour MRINV (voir MathLibC.h)                    //
//                                                                          //
// A : Matrice � inverser repr�sentant le syst�me d'equations lin�aires     //
// H : Pas calcul� par r�solution du syst�me d'�quations lin�aires          //
// Y : Vecteur solution mis � jour                                          //
// Evol : Vecteur repr�sentant l'�volution de la solution                   //
//        d'une it�ration sur l'autre                                       //
//                                                                          //
// F_relax : Valeur du syst�me au cours de la relaxation                    //
//                                                                          //
/*--------------------------------------------------------------------------*/
void MethodesNumeriques::newton_raph(int itemax, int* ite, double* crit_conver, double* crit_arret,
	double* sol, int derivee, double pas, int relax, enum STATUT* statut, ModeleEANL* modele, double* Fsol)
{
	int i, ite_relax;
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

	// Definition du probleme
	int ndim = modele->Get_dimension();

	H = new double[ndim + 1];
	Y = new double[ndim + 1 ];
	Evol = new double[ndim + 1];
	F_relax = new double[ndim + 1 ];

	A = new double*[ndim+1];
	for (i = 1; i <= ndim; i++)
	{
		A[i] = new double[ndim + 1 + 1];
	}


	// D�but de la boucle de calcul - Initialisation
	*ite = 0;
	modele->EvaluerF(sol,Fsol); // Intitialisation de F
	Norme = MethodesNumeriques::norme(Fsol, ndim); // Initialisation de la norme de F
	*statut = EN_COURS; // Initialisation du statut de la m�thode

	do
	{
		//  Test sur le crit�re de convergence
		if (Norme > *crit_conver)
		{
			(*ite)++; // Incr�ment du nombre d'it�rations
			if (*ite <= itemax)
			{
				// Calcul du Jacobien
				switch (derivee)
				{
				case 1: // Jacobien analytique
					modele->EvaluerJacob(sol,A);

					break;
				case 2: // Jacobien num�rique
					MethodesNumeriques::JAC_NUM(sol, Fsol, A, modele,pas);
					break;
				}
				// Remplissage de la matrice A : d�finition du syst�me d'�quations lin�aires
				for (i =1; i <= ndim; i++)
				{
					A[i][ndim+1] = -Fsol[i];
				}

				// Appel a Mrinv : proc�dure de r�solution de systeme d'equations lineaires
				MethodesNumeriques::MRINV(A, A, ndim, 0, Deter, pivot_min, H, type_calcul);

				// Calcul de sol au pas suivant
				switch (relax)
				{
				case 1: // Pas de ralaxation
					for (i = 1; i <= ndim; i++)
					{
						Y[i] = sol[i] + H[i];
					}
					break;
				case 2: // Relaxation num�rique
					alpha = 1; // Initialisation du coefficient alpha
					ite_relax = 0;
					do {
						for (i =1; i <= ndim; i++)
						{
							Y[i] = sol[i] + alpha * H[i]; // Calcul de la nouvelle solution
						}
						modele->EvaluerF(Y, F_relax); // Evaluation de F
						alpha = alpha / 2.; // Relaxation par r�duction du param�tre alpha
						ite_relax++; // On incr�mente le compteur
					} while (norme(F_relax, ndim) > Norme && ite_relax < 5); // Tant que la norme du nouveau F (F_relax)
																			 //est sup�rieure � l'ancienne

				}


				// Calcul de l'�volution entre deux it�rations
				for (i = 1; i <= ndim; i++)
				{
					Evol[i] = fabs((Y[i] - sol[i]) / sol[i]);
				}

				// Si la m�thode ne stagne pas, on met � jour la solution
				if (maxi(Evol, ndim) > *crit_arret)
				{
					for (i = 1; i <= ndim; i++)
					{
						sol[i] = Y[i];
					}
					modele->EvaluerF(sol, Fsol); // Mise � jour de F
					Norme = norme(Fsol, ndim); // Mise � jour de la norme de F
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

	} while (*statut == EN_COURS);

	// Mise � jour des param�tres de sortie
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



/*--------------------------------------------------------------------------*/
/* D�finition du jacobien num�rique										*/
/*--------------------------------------------------------------------------*/
void MethodesNumeriques::JAC_NUM(double* x,double* f,double**jacob, ModeleEANL* modele,double pas)
{
	int i, j;
	double* g;
	double u;

	int ndim = modele->Get_dimension();

	// Allocation dynamique de g
	// -------------------------
	g = new double[ndim + 1];

	for (j = 1; j <= ndim; j++)
	{
		u = x[j];
		x[j] = x[j] + x[j] * pas;

		for (i = 1; i <= ndim; i++)
		{
			modele->EvaluerF(x,g);
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
double MethodesNumeriques::norme(double *x, int ndim)
{
	int i;
	double norm = 0;

	for (i = 1; i <= ndim; i++)
	{
		norm = norm + pow(x[i], 2);
	}

	norm = sqrt(norm);
	return (norm);
}

/*----------------------------------------------------------------------------*/
/* Fonction visant � rechercher la composante maximal d'un vecteur de r�els  */
/*----------------------------------------------------------------------------*/
double MethodesNumeriques::maxi(double* Tab, int ndim)
{
	int i;
	double max;

	max = Tab[1];

	for (i = 1; i <= ndim; i++)
	{
		if (Tab[i] > max)
		{
			max = Tab[i];
		}
	}
	return (max);
}



void MethodesNumeriques::MRINV(double **A, double **B, int N, int NRC,
	double &DETER, double EPS,double *X, int INDIC)

	/*---------------------------------------------------------------------
	  R�le : M�thode de GAUSS JORDAN avec PIVOT MAXIMUM pour la r�solution
			 d'un syst�me de N �quations lin�aires ou l'inversion d'une
			 matrice
	  ---------------------------------------------------------------------
	  Arguments en entree :

		  - A     : matrice des coefficients augment�e du deuxi�me membre
					dans la N+1 �me colonne
		  - N     : nombre d'�quations ou dimension de la matrice B inverser
		  - NRC   : dimension de B et B (>= N)
		  - EPS   : plus petite valeur acceptable pour un pivot (en valeur
					absolue)
		  - INDIC : NEGATIF, calcul de la matrice inverse de B
					NUL, calcul de la solution du syst�me et de l'inverse de
					la matrice des coefficients
					POSITIF, r�solution du syst�me seulement

	  Arguments en sortie

		  - B     : Matrice contenant la matrice inverse apr�s traitement
					*** si on ne souhaite pas conserver B, appeler le SP par :
							MRINV(B,B,N,NRC,DETER,EPS,X,INDIC)
					La solution est calcul�e dans la N+1�me colonne de B puis
					rang�e dans X,
		  - X     : vecteur solution,
		  - DETER : valeur du determinant de la matrice des coeffcients

	  =======================================================================*/
{

	int I, J, K;
	int MAX, KM1, ISCAN, JSCAN;
	int IROWK, JCOLK, IROWI, JCOLI, INTCH, NM1, IP1;
	int JTEMP, IROWJ, JCOLJ;
	double PIVOT, AIJCK;
	int DUMMY;

	double *Y;     // (Dim 100)
	int *IROW;
	int *JCOL;
	int *JORD;

	/* allocation des matrices de travail
	   ---------------------------------- */
	Y = new double[N + 1];
	IROW = new int[N + 1];
	JCOL = new int[N + 1];
	JORD = new int[N + 1];

	MAX = N;
	if (INDIC >= 0) MAX = N + 1;
	//      if (N>100) goto E990 ;
	for (I = 1; I <= N; I++)
		for (J = 1; J <= MAX; J++)
			B[I][J] = A[I][J];
	DETER = 1.0;

	/* DEBUT DE LA PROCEDURE D ELIMINATION
	   ----------------------------------- */
	for (K = 1; K <= N; K++)
	{
		KM1 = K - 1;
		PIVOT = 0.0;
		for (I = 1; I <= N; I++)
		{
			for (J = 1; J <= N; J++)
			{
				if (K == 1) goto E9;
				for (ISCAN = 1; ISCAN <= KM1; ISCAN++)
				{
					for (JSCAN = 1; JSCAN <= KM1; JSCAN++)
					{
						if (I == IROW[ISCAN]) goto E11;
						if (J == JCOL[JSCAN]) goto E11;
					}
				}
			E9:             if (fabs(B[I][J]) <= fabs(PIVOT)) goto E11;
				PIVOT = B[I][J];
				IROW[K] = I;
				JCOL[K] = J;
			E11:			DUMMY = 0; // CONTINUE
			}
		}

		if (fabs(PIVOT) < EPS) goto E980;
		IROWK = IROW[K];
		JCOLK = JCOL[K];
		DETER = DETER * PIVOT;
		for (J = 1; J <= MAX; J++)
		{
			B[IROWK][J] = B[IROWK][J] / PIVOT;
		}
		B[IROWK][JCOLK] = 1.0 / PIVOT;
		for (I = 1; I <= N; I++)
		{
			AIJCK = B[I][JCOLK];
			if (I == IROWK) goto E18;
			B[I][JCOLK] = -AIJCK / PIVOT;
			for (J = 1; J <= MAX; J++)
			{
				if (J != JCOLK) B[I][J] = B[I][J] - AIJCK * B[IROWK][J];
			}
		E18:		 DUMMY = 0;  // CONTINUE		  
		}
	}

	/* ORDONNER LE VECTEUR SOLUTION
	   ---------------------------- */
	for (I = 1; I <= N; I++)
	{
		IROWI = IROW[I];
		JCOLI = JCOL[I];
		JORD[IROWI] = JCOLI;
		if (INDIC >= 0) X[JCOLI] = B[IROWI][MAX];
	}

	/* SIGNE DU DETERMINANT
	   -------------------- */
	INTCH = 0;
	NM1 = N - 1;
	for (I = 1; I <= NM1; I++)
	{
		IP1 = I + 1;
		for (J = IP1; J <= N; J++)
		{
			if (JORD[J] >= JORD[I]) goto E22;
			JTEMP = JORD[J];
			JORD[J] = JORD[I];
			JORD[I] = JTEMP;
			INTCH = INTCH + 1;
		E22:		DUMMY = 0;  // CONTINUE
		}
	}
	if (INTCH / 2 * 2 != INTCH) DETER = -DETER;

	/* REMISE EN ORDRE DE LA MATRICE INVERSE
	   ------------------------------------- */
	if (INDIC > 0) goto E900;
	for (J = 1; J <= N; J++)
	{
		for (I = 1; I <= N; I++)
		{
			IROWI = IROW[I];
			JCOLI = JCOL[I];
			Y[JCOLI] = B[IROWI][J];
		}
		for (I = 1;I <= N; I++)
		{
			B[I][J] = Y[I];
		}
	}
	for (I = 1; I <= N; I++)
	{
		for (J = 1; J <= N; J++)
		{
			IROWJ = IROW[J];
			JCOLJ = JCOL[J];
			Y[IROWJ] = B[I][JCOLJ];
		}
		for (J = 1; J <= N; J++)
		{
			B[I][J] = Y[J];
		}
	}
E900: goto E1000;

E980: ErreurMathLibC(6);

E1000:delete[] Y;
	delete[] IROW;
	delete[] JCOL;
	delete[] JORD;  // END ;
}


void MethodesNumeriques::ErreurMathLibC(int Numero)
{
	switch (Numero)
	{
	case  1: cout << "Erreur TRIDIA 01 : Le vecteur A n'est pas alloue !\n";
		break;
	case  2: cout << "Erreur TRIDIA 02 : Le vecteur B n'est pas alloue !\n";
		break;
	case  3: cout << "Erreur TRIDIA 03 : Le vecteur C n'est pas alloue !\n";
		break;
	case  4: cout << "Erreur TRIDIA 04 : Le vecteur D n'est pas alloue !\n";
		break;
	case  5: cout << "Erreur TRIDIA 05 : La dimension N est inf�rieure � 2 !\n";
		break;
	case  6: cout << "Erreur MRINV 01 : Pivot nul -> Matrice singuli�re !\n";
		break;
	}
	cout << " -> Programme avorte \n";
	exit(0);
}