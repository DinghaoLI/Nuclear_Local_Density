/**
 * @file Density.cpp
 *
 * Fichier de la calcul pour le calcul de la densité
 */

#include "Density.h"

/** Constructeur
 *
 */
Density::Density()
{
}

/** Constructeur avec des valeurs données
 *
 * @param _basis
 *
 */
Density::Density(Basis _basis)
{
    basis = _basis;
}


/**
 *
 * Cette fonction permet de calculer la densité
 *
 * @param rVals Les valeurs sur l’axe R
 * @param zVals Les valeurs sur l’axe Z
 *
 * @return la matrice de la densité
 */
arma::mat Density::calcDensity(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = 0;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            arma::mat funcA = basis.basisFunc(m,  n,  n_z, zVals, rVals);
                            arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
                            result += funcA % funcB * rho(ia,ib); // mat += mat % mat * double
                            // std::cout << "Basis vector " << ia << " " << ib << ": m=" << mp << " n=" << np << " n_z=" << n_zp << std::endl;
                            ib++;
                        }
                    }
                }
                ib = 0;
                ia++;
            }
        }
    }
    return result;
}


/**
 *
 * Cette fonction permet de calculer la densité de façon très rapide mais il y a un ecart avec l'algorithme direct
 *
 * @param rVals Les valeurs sur l’axe R
 * @param zVals Les valeurs sur l’axe Z
 *
 * @return la matrice de la densité
 */
arma::mat Density::calcDensityOp(arma::vec rVals, arma::vec zVals)
{

    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = -1;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    arma::mat sum = arma::zeros(rVals.n_rows, zVals.n_rows);
    arma::imat sumN_zMax = arma::sum(basis.n_zMax, 1);
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            arma::vec rPartA = basis.rPart(rVals, m, n);
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::mat funcA = rPartA * basis.zPart(zVals, n_z).t();
                sum = arma::zeros(rVals.n_rows, zVals.n_rows);
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    if (m != mp)
                    {
                        ib += sumN_zMax(mp);
                        continue;
                    }
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        arma::vec rPartB = basis.rPart(rVals, mp, np);
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            ib++;
                            arma::mat funcB = rPartB * basis.zPart(zVals, n_zp).t();
                            result += funcA % funcB * rho(ia,ib); // mat += mat % mat * double
                        }
                    }
                }
                ib = -1;
                ia++;
            }
        }
    }
    return result;
}

/**
 *
 * Cette fonction permet de calculer la densité, Il est assez vite sans biais
 *
 * @param rVals Les valeurs sur l’axe R
 * @param zVals Les valeurs sur l’axe Z
 *
 * @return la matrice de la densité
 */
arma::mat Density::calcDensity1(arma::vec rVals, arma::vec zVals)
{
    arma::mat rho;
    rho.load("rho.arma", arma::arma_ascii);

    int ia = 0;
    int ib = -1;
    arma::mat result = arma::zeros(rVals.n_rows, zVals.n_rows); // number of points on r- and z- axes
    arma::mat sum = arma::zeros(rVals.n_rows, zVals.n_rows);
    arma::imat sumN_zMax = arma::sum(basis.n_zMax, 1);
    for (int m = 0; m < basis.mMax; m++)
    {
        for (int n = 0; n < basis.nMax(m); n++)
        {
            arma::vec rPartA = basis.rPart(rVals, m, n);
            for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
            {
                arma::mat funcA = rPartA * basis.zPart(zVals, n_z).t();
                sum = arma::zeros(rVals.n_rows, zVals.n_rows);
                for (int mp = 0; mp < basis.mMax; mp++)
                {
                    if (m != mp)
                    {
                        ib += sumN_zMax(mp);
                        continue;
                    }
                    for (int np = 0; np < basis.nMax(mp); np++)
                    {
                        arma::vec rPartB = basis.rPart(rVals, mp, np);
                        for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++)
                        {
                            ib++;
                            if (n_z > n_zp)
                            {
                                arma::mat funcB = rPartB * basis.zPart(zVals, n_zp).t();
                                sum += 2 * funcB * rho(ia,ib); // mat += mat % mat * double
                            }
                            else if (n_z == n_zp)
                            {
                                arma::mat funcB = rPartB * basis.zPart(zVals, n_zp).t();
                                sum += funcB * rho(ia,ib); // mat += mat % mat * double
                            }
                        }
                    }
                }
                result += funcA % sum;
                ib = -1;
                ia++;
            }
        }
    }
    return result;
}