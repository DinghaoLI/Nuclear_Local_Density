#ifndef DEF_POLY
#define DEF_POLY
/**
 * @file Poly.h
 *
 * Interface de la classe Poly
 */
#include <iostream>
#include <armadillo>

/**
 * @class Poly
 *
 * Permet de calculer le polynôme d'Hermite et Laguerre généralisé
 *
 * ## Polynôme d'Hermite
\f[
\begin{eqnarray*}
H_0(\zeta)&=&1\\
H_1(\zeta)&=&2\zeta\\
H_{n}(\zeta)&=&2\zeta H_{n-1}(\zeta)-2(n-1)H_{n-2}(\zeta)
\end{eqnarray*}
\f]

## Polynôme de Laguerre généralisé

\f[
\begin{eqnarray*}
L^{m}_0(\eta)&=&1\\
L^{m}_1(\eta)&=&1+m-\eta\\
L^{m}_{n}(\eta)&=&\left(2+\frac{m-1-\eta}{n}\right)L^{m}_{n-1}(\eta)-\left(1+\frac{m-1}{n}\right)L^{m}_{n-2}(\eta)
\end{eqnarray*}
\f]
 */
class Poly
{
private:
    arma::mat polyHermite;
    arma::cube polyLaguerre;
public:
    Poly(void);
    void calcHermite(int, arma::vec);
    arma::vec hermite(int);
    void calcLaguerre(int, int, arma::vec);
    arma::vec laguerre(int, int);
};

#endif