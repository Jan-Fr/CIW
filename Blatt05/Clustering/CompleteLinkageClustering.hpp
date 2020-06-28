#pragma once


#include <iostream>
#include <stdlib.h>


#include "Naomini/Forward.hpp"
#include "Naomini/Helpers.hpp"



namespace Naomini {
double** MakeDistanceMatrix(MoleculeVector mols);

std::tuple<unsigned int, unsigned int> GetIndices(double** M, double sim, double greaterThan, unsigned int n);


double** OverwriteMatrix(double** M, unsigned int low_i, unsigned int low_j, unsigned int n);

void AppendToCluster(std::vector<std::vector<unsigned int>> &cluster, unsigned int low_i, unsigned int low_j);

} // end namespace Naomini

