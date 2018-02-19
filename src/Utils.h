#ifndef DEF_UTILS
#define DEF_UTILS
/**
 * @file Utils.h
 *
 * Interface de la classe Utils
 *
 */
#include <iostream>
#include <armadillo>

/**
 * @class Utils
 *
 * les fonctions de calculation mathématique
 */
class Utils
{
public:
    static arma::mat derivative(arma::mat, arma::mat, arma::mat, arma::mat);
    static int fac(int);
    static void matToFile(arma::mat, std::string);
    static std::string cubeToRaw(const arma::cube&);
};

#endif