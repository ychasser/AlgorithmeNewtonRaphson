/* ================================================================ *
 *  Fichier:     MathLibC.h					    *					
 *  Description: DLL de fonctions mathematiques pour les TP AFP	    *
 *               de 1ere et 2eme Annee GI                           *							
 *  Rev:         1.0				                    *								
 *  Date:        20 Novembre 2000				    *					
 *  Auteur:      Gilles HETREUX					    *					
 *  Mail:        Gilles.Hetreux@ensigct.fr			    *																	
 * ---------------------------------------------------------------- *
 *  INPT/ENSIGC et LGC(UMR CNRS 5503), Equipe AFP et GI		    *
 *  18, chemin de la Loge				            *					
 *  31058 TOULOUSE cedex, FRANCE				    *				
 *  Tel (33) 05.62.25.23.95					    *																
 * ================================================================ */



/*=====================================================================*/
  void __declspec(dllexport) TRIDIA ( double *A, double *B, double *C,
				      double *D, int N ) ;
/*---------------------------------------------------------------------
  Rôle : Ce sous-programme résoud le système linéaire suivant :
      B(1)*X(1)+C(1)*X(2)                                        =D(1)
      A(2)*X(1)+B(2)*X(2)+C(2)*X(3)                              =D(2)
                A(3)*X(2)+B(3)*X(3)+C(3)*X(4)                    =D(3)
                 ..........................
                                            A(N)*X(N-1)+B(N)*X(N)=D(N)
  ---------------------------------------------------------------------
  Arguments en entree :
	- A,B,C : vecteurs contenant les éléments tridiagonaux,
	          -> indicés de 1 à N. 
	- D	: vecteurs contenant le second membre,
	          -> indicé de 1 à N.
	- N	: nombre de ligne (ou colonne) de la matrice tridiagonale

  Arguments en sortie :
    	- A,B,C : non modifiés,
	- D	: vecteur solution du système lineaire
  =======================================================================*/
  
  
  
  
/*=====================================================================*/
 
   void __declspec(dllexport) NBOR (double Xmin, double Xmax, 
                                    double Precision, double &Critere, 
				    double &Xsol, double (*Fonction)(double) ) ;
/*---------------------------------------------------------------------
  Rôle : Ce programme recherche le minimum d'une fonction unimodale 
         d'une variable suivant l'algorithme du nombre d'or   
  ---------------------------------------------------------------------
  Arguments en entree :
	- Xmin     : borne inféreure du domaine de recherche
	- Xmax	   : borne supérieure du domaine de recherche
	- Précision: facteur de reduction de l'intervalle initiale
	- Fonction : fonction à minimiser. Elle doit etre unimodale
	             (1 seul optimum) sur l'intervalle

  Arguments en sortie
	- Critere : valeur du critere a l'optimum
	- Xsol : valeur de l'optimum
  =======================================================================*/



/*=====================================================================*/

  void __declspec(dllexport) MRINV ( double **A, double **B, int N, int NRC, 
				     double &DETER, double EPS, double *X, 
				     int INDIC ) ;

/*---------------------------------------------------------------------
  Rôle : Méthode de GAUSS JORDAN avec PIVOT MAXIMUM pour la résolution
         d'un système de N équations linéaires ou l'inversion d'une 
	 matrice
  ---------------------------------------------------------------------
  Arguments en entree :
      - A     : matrice des coefficients augmentée du deuxième membre
                dans la N+1 ème colonne
      - N     : nombre d'équations ou dimension de la matrice a inverser
      - NRC   : dimension de A et B (>= N)
      - EPS   : plus petite valeur acceptable pour un pivot (en valeur
                absolue)
      - INDIC : NEGATIF, calcul de la matrice inverse de A
                NUL, calcul de la solution du système et de l'inverse de 
                la matrice des coefficients
                POSITIF, résolution du système seulement

  Arguments en sortie
      - B     : Matrice contenant la matrice inverse aprés traitement 
                *** si on ne souhaite pas conserver A, appeler le SP par :
                         MRINV(A,A,N,NRC,DETER,EPS,X,INDIC)
                La solution est calculée dans la N+1ème colonne de B puis 
                rangée dans X,
      - X     : vecteur solution,
      - DETER : valeur du determinant de la matrice des coeffcients
  =======================================================================*/





/*=====================================================================*/

  void __declspec(dllexport) IDENMU(double *XOB,double *XMO,int &N,double *AI,
				    int &NPA,double **TAB,double **HESS,
				    double &OM,double &ECAR, int &KIMP,
				    int &NAP,double &ZMIN, 
				    void (* MODELE) (double *, double *) ) ;

/*--------------------------------------------------------------------- 
   Role :
     CE SOUS-PROGRAMME UTILISE LA METHODE DE GAUSS-NEWTON POUR
     IDENTifIER LES VALEURS DES PARAMETRES D'UN MODELE QUI MINIMISENT
     LE CRITERE D'ECART: SOMME DE (XOB(I)-XMO(I,AI))**2
     XOB(I) SONT DES DONNEES ET XMO(I,AI) SONT LES VALEURS  FOURNIES
     PAR LE MODELE POUR UN ENSEMBLE AI DE VALEURS DES PARAMETRES.
     LES VALEURS DE XMO(I,AI) SONT OBTENUES PAR APPEL DU SOUS-PROGRAMME
     MODELE:  CALL MODELE(AI,XMO).

     XOB(N)        VECTEUR DES DONNEES
     XMO(N)        VECTEUR DE TRAVAIL.AU RETOUR IL CONTIENT LES
                   VALEURS CALCULEES PAR LE MODELE AU MINIMUM DU
                   CRITERE.
     N             NOMBRE DE POINTS DE DONNEES.
     AI(NPA)       VECTEUR CONTENANT
                   A L'APPEL:ESTIMATION INITIALE DES PARAMETRES
                   AU RETOUR:VALEURS DE CEUX-CI AU MINIMUM.
     NPA           NOMBRE DE PARAMETRES A IDENTifIER.
     TAB(NPA,N)    TABLEAU DE TRAVAIL.
     HESS(NPA,NPA):TABLEAU DE TRAVAIL.AU RETOUR IL CONTIENT LA
                   MATRICE DE LA FORME QUADRATIQUE ASSOCIEE AU CRITERE
                   AU MINIMUM.
     OM            FACTEUR DE RELAXATION (VALEUR POSSIBLE:0.5).
     ECAR          VARIATION RELATIVE DU CRITERE ENTRE DEUX ETAPES
                   EN DESSOUS DE LAQUELLE LA RECHERCHE EST ARRETEE.
     KIMP          IMPRESSION DES RESULTATS A CHAQUE PAS SI KIMP=1.
     NAP           NOMBRE D'APPELS DU MODELE PENDANT LA RECHERCHE.
     ZMIN          VALEUR DU CRITERE AU MINIMUM.

 =======================================================================*/





/********************************************************************/

   void __declspec(dllexport)ROSEN (int NO,int N, int M,int IPRINT, double EPS,
                                    double EPS1, double *X,double *YX, 
				    double **GX, double *G, int &NAF, 
				    double *ACC, int ICHOIX, int ITMAX,
				    int IFOIS, 
				    double (* FCT) (double *, int, int ),
				    void (* GRADF) (double*, double *, int),
				    void (* CONS)(double *,double *,double **,int,int),
				    double &FMAX, double ECD, double XINF) ;

/***********************************************************************
C
C  Role :
C     Maximise une fonction non lineaire en presence de contraintes
C     lineaires du type G(X)<=0 ou G(X)=0
C
C     Attention: 1-Pour minimiser une fonction,il faut changer
C                  signe.
C                2-Le programme est ecrit en double precision.
C                3-Dans cette version,le nombre de variables et
C                  de contraintes sont limites a 50.Si l'on desire
C                  augmenter cette valeur,il faut changer les
C                  dimensions des tableaux de travail:
C                   - vecteurs ICC,BB,DIR,GG,GGG,G0,ZERO dans ROSEN
C                   - vecteur   AA dans POINT
C                   - vecteur   W  dans MINOR2
C
C     La procedure utilise la methode du gradient projete de Rosen
C     avec recherche du pas optimal par la methode du nombre d'or.
C
C     Entrees
C     -------
C
C     NO           No  du fichier de sortie(par defaut NO=6).
C     N            Nombre de variables independantes (<=50 ).
C     M            Nombre de contraintes(<=50 ).
C     IPRINT      =0 pas d'impression.
C                 =1 impressions intermediaires.
C     EPS          Precision en dessous de laquelle une contrainte
C                  est  consideree saturee(par defaut 1.e-3).
C     EPS1         Provoque l'arret de la recherche lorsque la norme de
C                  la direction de montee est inferieure a EPS1
C                  (Valeur par defaut 1.e-4).
C     X(N)         Valeurs initiales des variables,qui doivent verifier
C                  les contraintes.
C     YX(N)        Gradient du critere (a definir dans la subroutine GRAD).
C     GX(M,N)      Jacobien des contraintes(a definir dans la subroutine CONS.
C     G(M)         Valeurs des contraintes(a definir dans CONS).
C     ECD          Amplitude de l'intervalle contenant la solution lors
C                  de la recherche de la longueur optimale du pas dans le
C                  sous-programme MINOR2(recherche monodimensionnelle)
C                  (valeur par defaut 1.e-4).
C     XINF         Borne maximale admissible pour la recherche du pas
C                  (valeur par defaut 10).
C     ITMAX        Nombre maximal d'iterations (valeur par defaut 50).
C     ifOIS        Si IPRINT=1,il y a des impressions intermediaires
C                  toutes les ifOIS iterations (valeur par defaut 1).
C     ICHOIX       Mode de calcul du gradient
C                  <0 : gradient analytique
C                  =0 : dif ferences finies du premier ordre du type
C                       X(I)+ACC(I)
C                  >0 : dif ferences finies du second ordre du type
C                       X(I)+ACC(I) puis X(I)-ACC(I).
C     ACC(N)       Vecteur des accroissements pour le calcul du gradient
C                  par differences finies.
C                  -------------------
C     FCT          Nom du sous-programme de type function dans lequel
C                  l'utilisateur definit la fonction a maximiser.
C                  A declarer EXTERNAL.
C                  -------------------
C     GRADF        Nom du sous-programme  de type subroutine dans lequel
C                  l'utilisateur definit le gradient analytique du critere
C                  (ce qui correspond a ICHOIX<0).Dans ce cas on peut
C                  initialiser le vecteur ACC avec des valeurs quelconques.
C                  Dans le cas contraire (ICHOIX>=0) il faut initialiser
C                  ICHOIX avec des valeurs assez faibles (0.05 pour tout
C                  i=1 a n par exemple), et on doit ecrire un sous programme
C                  GRADF vide :
C                  SUBROUTINE GRADF(X,YX,N)
C                  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                  DIMENSION X(N),YX(N)
C                  RETURN
C                  END
C                  Il faut declarer GRADF EXTERNAL
C                  -------------------
C     CONS         Nom du sous-programme dans lequel l'utilisateur defi-
C                  nit les contraintes et eventuellement leurs derivees
C                  A declarer EXTERNAL
C                  --------------------
C
C     Sorties
C     -------
C
C     X(N)         Valeurs des variables en fin de recherche(point
C                  solution)
C     FMAX         Valeur du critere en fin de recherche.
C     NITER        Nombre d'iterations.
C     NAF          Compteur laisse a la disposition de l'utilisateur
C                  (voir ci-dessous) .
C
C     Variables et tableaux de travail
C     --------------------------------
C
C     ICC(M)       Indices des contraintes saturees(surdim. 50).
C     AA(M*M)      Coefficients du systeme lineaire donnant les
C                  parametres de Lagrange correspondant aux contraintes
C                  saturees(surdim. 2500).
C     BB(M)        Valeurs des parametres de Lagrange(surdim. 50).
C     DIR(N)       Direction de montee(surdim. 50).
C     GG,GGG,W     Tableaux de stockage(surdim. 50).
C     ZERO,G0      Tableaux de travail(surdim. 50).
C
C     Sous-programmes appeles
C
C     -----------------------------------------------------------------
C     POINT,DCHOLE,MINOR2,GRADFC : Sous-programmes de la bibliotheque
C                                  standard.
C     -----------------------------------------------------------------
C     FCT(X,NAF,N)        Sous-programme de type FUNCTION ou l'utilisa-
C                         teur definit le critere.
C                         NAF est laisse a la disposition de l'utilisateur
C                         qui peut par exemple l'utiliser pour compter le
C                         nombre d'appels au critere au cours de la
C                         recherche (en ecrivant dans la fonction
C                         NAF=NAF+1).
C                         A declarer EXTERNAL dans le programm'appelant.
C     GRADF(X,YX,N)        Sous-programme de type SUBROUTINE ou l'utilisa-
C                         teur definit le gradient analytique du critere
C                         (voir la remarque ci-dessus).
C                         A declarer EXTERNAL dans le programme appelant.
C     CONS(X,G,GX,N,M)    Sous-programme  de type SUBROUTINE ou l'utilisateur
C                         definit les contraintes et eventuellement leurs
C                         derivees.Les contraintes doivent toutes etre
C                         ecrites sous forme <=0.Par ailleurs pour reduire
C                         le temps de calcul,il est preferable de definir
C                         les derivees des contraintes (qui sont constantes)
C                         dans le programme appelant.
C                         A declarer EXTERNAL dans le programme appelant.
C
C     Remarque            Dans la mesure ou ces trois sous-programmes
C     --------            sont declares de type EXTERNAL,l'utilisateur
C                         peut changer leurs noms.
C
C     -----------------------------------------------------------------
C
C
C********************************************************************/



/* -------------------------------------------------------------------- */

     void GRADFC ( double *X, double *YX,int N,int ICHOIX, double *ACC, 
                   int &NAF, void (* GRADF) (double *, double *, int),
		   double (* FCT) (double *, int , int ) ) ;

/* --------------------------------------------------------------------
C
C     DEFINITION DU GRADIENT DU CRITERE
C
C     ENTREES
C     -------
C     X       (VOIR MODULE ROSEN)
C     N       (VOIR MODULE ROSEN)
C     ICHOIX  (VOIR MODULE ROSEN)
C     ACC     (VOIR MODULE ROSEN)
C     NAF     (VOIR MODULE ROSEN)
C     GRADF   (VOIR MODULE ROSEN)
C     FCT     (VOIR MODULE ROSEN)
C
C     SORTIE
C     ------
C     YX      VECTEUR GRADIENT DU CRITERE
C
C     SOUS-PROGRAMMES APPELES
C     -----------------------
C
C     SUBROUTINE UTILISATEUR GRADF ET FONCTION UTILISATEUR FCT
C
C----------------------------------------------------------------------*/



/***********************************************************************/

   void POINTFUNC ( double *X,double *YX, double **GX,double *DIR, int *ICC,
		int N, int M, int L, int &JA, double EPS,double &C,
		int IPRINT, int NO, int IFOIS, int ICHOIX, int &NITER,
		double *ACC, int &NAF, 
		void (* GRADF) (double *, double *, int ),
		double (* FCT) (double *, int, int ) ) ;

/***********************************************************************
C
C
C     Ce module calcule le gradient projete et effectue un test
C     d'arret lorsque sa norme est arbitrairement faible
C
C     entree :
C     -------
C     X       (voir module Rosen)
C     GX      (voir module Rosen)
C     ICC     (voir module Rosen)
C     N       (voir module Rosen)
C     M       (voir module Rosen)
C     L       (voir module Rosen)
C     EPS     (voir module Rosen)
C     IPRINT  (voir module Rosen)
C     N0      (voir module Rosen)
C     IFOIS   (voir module Rosen)
C     ICHOIX  (voir module Rosen)
C     ACC     (voir module Rosen)
C     GRADF   (voir module Rosen)
C
C     sorties
C     -------
C     YX(N)   Gradient du critere (calcule par appel au module GRADFC)
C     DIR(N)  Gradient projete
C     JA      Indicateur d'arret . Si JA=0 on poursuit la recherche,
C             si JA=1 la recherche est arretee car on a atteint
C             une solution.
C     C       Inutilise dans cette version.
C
C     Sous-programmes appeles
C     -----------------------
C     Subroutines (bibliotheque) GRADFC et DCHOLE
C
C--------------------------------------------------------------*/




/************************************************************************/

    void MINOR2 ( double &ZMIN, double ZMAX,double ECD, double &CT,
		  double *X, double *DIR, int N,int NO,
		  double FCT (double *, int , int),
		  int IPRINT, int IFOIS,int &NITER,int &NAF) ;

/***********************************************************************
C
C     Recherche du pas optimal le long de la direction de
C     descente par la methode du nombre d'or
C
C     entrees
C     -------
C     ZMIN   Borne inferieure de l'intervalle de recherche
C     ZMAX   Borne superieure de l'intervalle de recherche
C     ECD    Precision de la recherche
C     X      (voir module ROSEN)
C     DIR    (voir module ROSEN)
C     N      (voir module ROSEN)
C     IPRINT (voir module ROSEN)
C     NO     (voir module ROSEN)
C     IFOIS  (voir module ROSEN)
C     NITER  (voir module ROSEN)
C     NAF    (voir module ROSEN)
C     FCT    (voir module ROSEN)
C
C     sorties
C     ------
C     ZMIN   Valeur optimale du pas
C     CT     Valeur optimale du critere (non utilisee dans cette version)
C
C----------------------------------------------------------------------*/

/******************************************************************/

    void DCHOLE ( double **A, int &M, double *B, int &N, int IOP,
	          int &IDEF, double &DET) ;

/******************************************************************/
