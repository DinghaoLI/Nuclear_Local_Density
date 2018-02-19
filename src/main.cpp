/**
 *@file main.cpp
 */

#include <fstream>
#include <armadillo>
#include "Poly.h"
#include "Basis.h"
#include "Density.h"
#include "time.h"

/**
 * @mainpage Calcul et tracé de la densité d'un système nucléaire
 *# Densité nucléaire locale
 *\f[\rho(\mathbf{r})\equiv \sum_a \sum_b \rho_{ab}\psi_a(\mathbf{r})\psi^*_b(\mathbf{r})\f]
 *avec \f[\rho_{ab}\f] une matrice donnée
 *et \f[\psi_a(\mathbf{r})\f] la fonction d'onde solution de l'équation de Schrödinger 3D
\f[
\psi_{m,n,n_z}(r_\perp, \theta, z)
        \equiv
    Z(z, n_z)
    .
    R(r_\perp, m, n)
    .
         e^{im\theta}
\f]
 * # Equation de Schrödinger
 * \f[\hat{H}_{(z)}\psi_n(z) = E_n\psi_n(z)\f]
 */

/**
 *La fonction principale
 *
 *Cette fonction a pour but de calculer les données nécessaires pour tracer les grraphes
 *
 *@return la fonction créer des fichier .txt et retourne 0
 */

int main()
{
    clock_t start, finish;
    double duration;

    Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);

    Density density(basis);

    std::cout<< "===Calcul de la fonction d'onde===" << std::endl;
    //Calcul de la densité 2D
    arma::vec rVals_psi = arma::linspace(-10,10,100);
    arma::vec zVals_psi = arma::linspace(-10,10,100);
    arma::mat psi = basis.basisFunc(0, 0, 0, zVals_psi, rVals_psi);
    Utils::matToFile(rVals_psi, "rVals_psi000.txt");
    Utils::matToFile(zVals_psi, "zVals_psi000.txt");
    Utils::matToFile(psi%psi, "Psi000.txt");

    psi = basis.basisFunc(0, 0, 1, zVals_psi, rVals_psi);
    Utils::matToFile(rVals_psi, "rVals_psi001.txt");
    Utils::matToFile(zVals_psi, "zVals_psi001.txt");
    Utils::matToFile(psi%psi, "Psi001.txt");

    psi = basis.basisFunc(0, 1, 1, zVals_psi, rVals_psi);
    Utils::matToFile(rVals_psi, "rVals_psi011.txt");
    Utils::matToFile(zVals_psi, "zVals_psi011.txt");
    Utils::matToFile(psi%psi, "Psi011.txt");

    psi = basis.basisFunc(1, 0, 1, zVals_psi, rVals_psi);
    Utils::matToFile(rVals_psi, "rVals_psi101.txt");
    Utils::matToFile(zVals_psi, "zVals_psi101.txt");
    Utils::matToFile(psi%psi, "Psi101.txt");
    std::cout<< "Psi.py pour le plot Psi" << std::endl;

    std::cout<< "===Calcul de la densité 2D===" << std::endl;
    //Calcul de la densité 2D
    arma::vec rVals_2D = arma::linspace(-10,10,100);
    arma::vec zVals_2D = arma::linspace(-10,10,100);
    arma::mat R = density.calcDensity1(rVals_2D, zVals_2D);
    Utils::matToFile(rVals_2D, "rVals.txt");
    Utils::matToFile(zVals_2D, "zVals.txt");
    Utils::matToFile(R, "plot2d.txt");
    std::cout<< "Density.py pour le plot 2D" << std::endl;



    std::cout<< "===Calcul de la densité 3D===" << std::endl;
    //Calcul de la densité 3D
    arma::vec zVals = arma::linspace(-20,20,64);

    arma::mat rVals;


    if (!rVals.load("rVals3d.txt"))
    {
        std::cout<< "=================================" << std::endl;
        std::cout<< "ERREUR: exécuter python rVals3d.py" << std::endl;
        std::cout<< "=================================" << std::endl;
    }

    start = clock();
    arma::cube results = arma::zeros(rVals.n_rows,rVals.n_cols,zVals.n_rows);
    arma::mat result = arma::zeros(rVals.n_rows,zVals.n_rows);

    for (uint i =0; i<rVals.n_cols; i++)
    {
        density.calcDensity1(rVals.col(i), zVals);
        for (uint k=0; k<results.n_slices; k++)
            results.slice(k).col(i) = result.col(k);
        finish = clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;
        std::cout<< i <<" :" << duration << "s" << std::endl;
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout<< duration << "s Temps de calcul de la densité 3D" << std::endl;

    //Enregistrement des données dans un fichier
    std::ofstream out("plot3d.raw");
    if (out.is_open())
    {
        out << Utils::cubeToRaw(results) << std::endl;
        out.close();
    }
    std::cout << "La matrice est stockée dans " << "plot3d.raw" << std::endl;

    //Analyse des performances des calculs
    std::cout<< "====Performances des calculs===" << std::endl;

    rVals = arma::linspace(-10,10,10);
    zVals = arma::linspace(-10,10,10);

    start = clock();

    arma::mat resOp = density.calcDensityOp(rVals, zVals);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout<< duration << "s après l'optimisation 1" << std::endl;

    start = clock();
    arma::mat res2 = density.calcDensity1(rVals, zVals);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout<< duration << "s après l'optimisation 2" << std::endl;

    start = clock();
    arma::mat res1 = density.calcDensity(rVals, zVals);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout<< duration << "s sans Optimization" << std::endl;

    std::cout << "Différence Op 1: " << arma::norm(res1 - resOp) << std::endl;
    std::cout << "Différence Op 2: " << arma::norm(res1 - res2) << std::endl;

    return 0;
}


