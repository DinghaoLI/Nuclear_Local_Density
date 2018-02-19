# Projet IPS-PROD - ENSIIE - 2018

Calcul et tracé de la densité d'un système nucléaire

## Pré-requis

* c++ 11
* cxxtest
* python
* armadillo
* doxygen
* valgrind

## Pour compiler le code et générer la documentation

```
make
```

## Pour formatter le code

```
make style
```

## Pour lancer les tests

```
make test
```

## Pour vérifier qu'il n'y a pas de fuites mémoires

```
make memcheck
```

## Pour supprimer les fichiers d'exécutions

```
make clean
```

## Pour la préparation de données (plot2d, plot3D, Psi)

Préparer les données pour dessiner les graphes.
```
python rVals3d.py
```
Puis, pour calculer les données, lancer la commande suivante:

```
./main
```

## Plot 2D

Librairies python requises: numpy, matplotlib, mpl_toolkits.mplot3d, scipy.io.
Après le calcul des données, lancer la commandes suivantes:

```
make plot2d
```

## Plot Psi

Sur votre python, Vous devez avoir les packages installés: numpy, matplotlib, mpl_toolkits.mplot3d, scipy.io.
Après la préparation de données, on a les données pour plot de Psi et on tape la commande suivante pour le montrer:

```
make plotpsi
```

## Plot 3D

Ouvrir avec blender le fichier plot3d.blend
Et cliquer sur F12 pour avoir la rendu en 3d.


