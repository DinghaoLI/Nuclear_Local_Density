/**
 * @file Utils.cpp
 *
 * Fichier d'implémentation pour la calculation mathématique
 */

#include "Utils.h"

/** Calcul de la factorielle de n
*
* Cette fonction a pour but de calculer la factorielle de n par récurrence.
*
* @param n entier
*
* @return le résultat de la factorielle de n
*/
int Utils::fac(int n)
{
    int f;

    if (n==0 || n==1)
        f=1;
    else
        f=fac(n-1)*n;

    return f;
}

/**
 * Calcul de dérivée
 *
 * Cette fonction a pour but d'approximer la dérivée du premier ordre
 *
 * @param F1 vecteur contenant les n valeurs d'une fonction
 * @param F2 vecteur contenant les n valeurs d'une fonction décalée d'une valeur
 * @param Z1 vecteur contenant les n points où sont évaluée la fonction
 * @param Z2 vecteur contenant les n points où sont évaluée la fonction décalée d'une valeur
 *
 * @return res vecteur contenant les valeurs de la fonction dérivée
 */
arma::mat Utils::derivative(arma::mat F1, arma::mat F2, arma::mat Z1, arma::mat Z2)
{
    arma::mat res = (F2-F1)/(Z2-Z1);
    return res;
}

/** Stockage des valeurs de la fonction d'onde
*
* Cette fonction a pour but de stocker les valeurs de la fonction d'onde dans un fichier txt
*
* @param Z
* @param title
*
* @return rien
*/
void Utils::matToFile(arma::mat Z, std::string title)
{
    std::ofstream out(title);
    if (out.is_open())
    {
        out << Z << std::endl;
        out.close();
    }
    std::cout << "La matrice est stockée dans " << title << std::endl;

}

/** La génération de données “Raw” pour dessiner le plot 3D
*
*
* @param &m la matrice cube contenant les données
*
* @return ss.str()
*/
std::string Utils::cubeToRaw(const arma::cube& m)
{
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);
    double theMin = 0.0;
    double theMax = m.max();
    for (uint k = 0; k < m.n_slices; k++)
    {
        for (uint j = 0; j < m.n_cols; j++)
        {
            for (uint i = 0; i < m.n_rows; i++)
            {
                uint v = 255 * (fabs(m(i, j, k)) - theMin) / (theMax - theMin);
                ss.put(v);
            }
        }
    }
    return ss.str();
}