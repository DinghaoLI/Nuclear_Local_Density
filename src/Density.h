#ifndef DEF_DENSITY
#define DEF_DENSITY
/**
 * @file Density.h
 *
 * Interface de la classe Density
 */
#include <iostream>
#include <armadillo>
#include <math.h>
#include <cmath>
#include "Poly.h"
#include "Utils.h"
#include "Basis.h"


/**
 * @class Density
 *
 * Permet de calculer la densité nucléaire locale
\f[\rho(\mathbf{r})\equiv \sum_a \sum_b \rho_{ab}\psi_a(\mathbf{r})\psi^*_b(\mathbf{r})\f]


avec \f[\rho_{ab}\f] une matrice donnée

et \f[\psi_a(\mathbf{r})\f] la fonction d'onde solution de l'équation de Schrödinger 3D
\f[
\psi_{a\equiv m,n,nz}(r, \theta, z)
        \equiv
    Z_a(z)
    .
    R_a(r)
    .
         e^{im\theta}
\f]

\f[\psi_a(\mathbf{r})\psi^*_b(\mathbf{r})\equiv \sum_a \sum_b \rho_{ab} R_a(r) Z_a(z) R_b(r) Z_b(z) \exp(im\theta - im\theta)\f]
\f[\psi_a(\mathbf{r})\psi^*_b(\mathbf{r})\equiv \sum_a \sum_b \rho_{ab} R_a(r) Z_a(z) R_b(r) Z_b(z)\f]
 */
class Density
{
private:
    Basis basis;

public:
    Density();
    Density(Basis);
    arma::mat calcDensity(arma::vec, arma::vec);
    arma::mat calcDensity1(arma::vec, arma::vec);
    arma::mat calcDensityOp(arma::vec, arma::vec);
};

#endif