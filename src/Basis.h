#ifndef DEF_BASIS
#define DEF_BASIS
/**
 * @file Basis.h
 *
 * Interface de la classe Basis
 */
#include <iostream>
#include <armadillo>
#include <math.h>
#include <cmath>
#include "Poly.h"
#include "Utils.h"

/**
 * @class Basis
 *
 * ## Troncature de la base
 * On a
 *\f[n_z^\textrm{max}(i) \equiv (N+2).Q^\frac{2}{3}+\frac{1}{2}-i.Q
\f]

et
\f[
m^\textrm{max} \equiv \textrm{sup}\left\{i:n_z^\textrm{max}(i)\ge 1\right\}.
\f]

Les valeurs des nombres quantique sont des **entiers** vérifiant
\f[\begin{eqnarray}
0 &\le m \le& m^\textrm{max}-1\\
0 &\le n \le& \frac{1}{2}(m^\textrm{max}-m-1)\\
0 &\le n_z \le& n_z^\textrm{max}(m+2n+1)-1.
\end{eqnarray}\f]

 *
 * Permet de calculer les bornes de la base
 */
class Basis
{
private:
    double br, bz;
    arma::vec zfactPart =
    {
        7.511255444649424828587030047762276930524e-1,
        5.311259660135984572385365242537567693773e-1,
        2.655629830067992286192682621268783846887e-1,
        1.08415633823009689789878853323749022452e-1,
        3.833071493144388678758035489700725876944e-2,
        1.212123635259875349071310579655202864027e-2,
        3.499099535541983943542762678345281363805e-3,
        9.351736874441379770202142628079336715494e-4,
        2.337934218610344942550535657019834178873e-4,
        5.510563799824823821067738395272340242989e-5,
        1.232199525075784976451417893101893685019e-5,
        2.627058214392003194960140445607747615225e-6,
        5.362460124872919207572643465049275268564e-7,
        1.051664954522700642243897278328369487056e-7,
        1.987459951592227331357006655358164809117e-8,
        3.628588825417294648844562498457222314549e-9,
        6.414499411475746158307919785557463603948e-10,
        1.100077573465546081464447454369726996765e-10,
        1.833462622442576802440745757282878327941e-11,
        2.97426912202769525933866781584586998244e-12,
        4.702732399958400033923267262424369697562e-13,
    };

    double calc_nzMax(int N, double Q, int i);
public:
    int mMax;
    arma::ivec nMax;
    arma::imat n_zMax;
    Basis();
    Basis(double, double, int, double);
    arma::vec zPart(arma::vec,int);
    arma::vec rPart(arma::vec, int, int);
    arma::mat basisFunc(int, int, int, arma::vec, arma::vec);
};

#endif