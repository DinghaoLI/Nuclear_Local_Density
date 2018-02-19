/**
 * @file Basis.cpp
 *
 * Fichier d'implémentation pour le calcul de la densité
 */

#include "Basis.h"


/** Constructeur des valeurs de la base
 *
 */
Basis::Basis()
{
}


/** Constructeur avec des valeurs données
 *
 * @param _br
 * @param _bz
 * @param N
 * @param Q
 */

Basis::Basis(double _br, double _bz, int N, double Q)
{
    bz = _bz;
    br = _br;

    int i = 0;
    int n, m;
    while (Basis::calc_nzMax(N,Q,i) >= 1)
    {
        i++;
    }
    mMax = i -1;
    arma::ivec vm = arma::linspace<arma::ivec>(0, mMax-1, mMax);
    nMax = (mMax - vm - 1) / 2 + 1;

    n_zMax.zeros(vm.n_rows, (mMax - 1) / 2 + 1);

    for (m=0; m < mMax; m++)
        for (n=0; n < nMax[m]; n++)
        {
            // n_zMax[m][n] =  (N + 2) * pow(Q, 2.0/3.0) + 0.5 - (m + 2*n) * Q;
            n_zMax(m, n) = Basis::calc_nzMax(N, Q, m + 2*n + 1);
        }


}

/** Calcul de la fonction d'nz_max
*
* On calcule nz_max par N, Q et i
* \f[n_z^\textrm{max}(i) \equiv (N+2).Q^\frac{2}{3}+\frac{1}{2}-i.Q\f]
*
* @param N
* @param Q
* @param i
*
* @return le résultat de calcul de nz_max
*/
double Basis::calc_nzMax(int N, double Q, int i)
{
    return (N + 2) * pow(Q, 2.0/3.0) + 0.5 - i * Q;
}



/** Calcul de la fonction de la partie sur Z
*
* On calcule par la formule suivante
* \f[ Z(z, n_z) \equiv \phi_{n_z}(z) = \frac{1}{\sqrt{b_z}} \frac{1}{\sqrt{2^{n_z} \sqrt{\pi}n_z!}} e^{-\frac{z^2}{2b_z^2}}H_{n_z}\left(\frac{z}{b_z}\right) \f]
*
* @param Z
* @param nz
*
* @return le résultat de la partie sur Z
*/
arma::vec Basis::zPart(arma::vec Z, int nz)
{
    Poly poly;
    double constPart;
    arma::vec expPart;

    constPart = (1 / std::sqrt(bz)) * zfactPart(nz);
    expPart =  arma::exp(-(Z % Z) / (2 * bz * bz));
    poly.calcHermite(nz, Z/bz);

    return constPart* expPart % poly.hermite(nz);
}


/** Calcul de la fonction de la partie sur R
*
* On calcule par la formule suivante
* \f[R(r_\perp, m, n) \equiv \frac{1}{b_{\perp}\sqrt{\pi}} \sqrt{\frac{n!}{(n+|m|)!}} e^{-\frac{r_{\perp}^2}{2b_{\perp}^2}} \left(\frac{r_{\perp}}{b_{\perp}}\right)^{|m|} L_n^{|m|}\left(\frac{r_{\perp}^2}{b_{\perp}^2}\right).\f]
*
* @param R
* @param m
* @param n
*
* @return le résultat de la partie sur R
*/
arma::vec Basis::rPart(arma::vec R, int m, int n)
{
    Poly poly;
    double constPart;
    arma::vec expPart;

    constPart = sqrt((double) Utils::fac(n)/ (double) Utils::fac(n+m)) / (br * sqrt(M_PI));
    expPart = arma::exp(-(R % R) / (2 * br * br));
    poly.calcLaguerre(m+1, n+1, R % R / (br * br));

    return constPart * expPart % pow(R / br, m) % poly.laguerre(m,n);
}



/** Calcul de la fonction de la Basis
*
* On calcule par la formule suivante
* \f[\psi_{m,n,n_z}(r_\perp, \theta, z) \equiv Z(z, n_z) . R(r_\perp, m, n) . e^{im\theta}\f]
*
* @param m
* @param n
* @param nz
* @param Z
* @param R
*
* @return le résultat de la Basis fonction
*/
arma::mat Basis::basisFunc(int m, int n, int nz, arma::vec Z, arma::vec R)
{
    return rPart(R, m, n) *  zPart(Z, nz).t();
}


