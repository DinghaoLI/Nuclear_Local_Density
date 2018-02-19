/**
 * @file Poly.cpp
 *
 * Fichier d'implémentation pour le calcul des polynômes d'Hermite et de Laguerre
 */

#include "Poly.h"

Poly::Poly(void)
{
}


/**
 * Calcul du polynôme d'Hermite
 *
 * Cette fonction a pour but de calculer le polynôme d'Hermite par récurrence avec des vecteurs
 *
 */
void Poly::calcHermite(int n, arma::vec Z)
{
    int i;

    if (n < 0)
    {
        std::cout << "Error !" << std::endl;
        exit(0);
    }

    polyHermite.ones(Z.n_rows,n+1);


    if (n == 0)
    {
        return;
    }

    polyHermite.col(1) = 2 * Z;

    for (i = 2; i <= n; i++)
    {
        polyHermite.col(i) = (2 * Z % polyHermite.col(i - 1)) - (2 * (i-1) * polyHermite.col(i - 2));
    }

}

arma::vec Poly::hermite(int n)
{
    return polyHermite.col(n);
}

/**
 * Calcul du polynôme de Laguerre par récurrence
 *
 * Cette fonction a pour but de calculer le polynôme de Laguerre par récurrence avec des vecteurs
 *
 * @param m
 * @param n
 * @param Z un vecteur de valeurs
 *
 * @return rien
 */
void Poly::calcLaguerre(int m,int n,arma::vec Z)
{
    int i;
    polyLaguerre.zeros(Z.n_rows,m+1,n+1);
    //cas n=0
    polyLaguerre.slice(0)++;
    //cas n=1
    arma::mat L1 = arma::zeros(Z.n_rows, m+1);
    arma::mat matZ = arma::zeros(size(L1));
    arma::mat M = arma::zeros(size(L1));
    for (i=0; i<=m; i++)
    {
        L1.col(i) = 1 + i - Z;
        matZ.col(i) = Z;
        M.col(i).fill(i);
    }
    polyLaguerre.slice(1) = L1;

    //les autres cas
    for (i=2; i<=n; i++)
    {
        polyLaguerre.slice(i) = (2 + (M - 1 - matZ) / i) % polyLaguerre.slice(i-1) - (1 + (M - 1)/ i) % polyLaguerre.slice(i-2);
    }

}

/**
 * Afficher le vecteur du polynôme de Laguerre
 *
 * Cette fonction a pour but d'obtenir le polynôme de Laguerre pour n,m donnés
 *
 * @param m
 * @param n
 *
 * @return un vecteur contenant les valeurs du polynôme de Laguerre pour n,m donnés
 */
arma::vec Poly::laguerre(int m,int n)
{
    if (polyLaguerre.n_cols < (unsigned) m || polyLaguerre.n_slices < (unsigned) n)
        throw "m,n out of Range";
    return polyLaguerre.slice(n).col(m);
}


