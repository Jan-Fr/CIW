#include <iostream>
#include <stdlib.h>

#include "Naomini/Forward.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeFactory.hpp"

#include "CompleteLinkageClustering.hpp"

using namespace Naomini;

/*Code compiliert, aber liefert falsche Ausgabe
Idee: 
1. Distanzmatrix erstellen
2. Kleinste Distanz finden und die Indices low_i und low_j bestimmen.
   Diese dann in einer Liste speichern.
3. Dabei wird keine neue kleinerer Distanzmatrix erstellt, sondern die
   alte wird ueberschrieben (und hat somit nach dem Ueberschreiben 2 identische Zeilen und Spalten)
4. Printen von zusammengehoerigen Clustern

Es fehlen noch:
Eine Methode, die doppelt vorkommende Cluster zusammenfuehrt
Eine Methode, die einzelne Cluster fuer Molekuele erstellt, deren Distanz
zu allen anderen Molekuelen den threshold ueberschreitet.
*/


int main(int argc, char* argv[]) {
  // check the number of arguments
  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " <molecule-file> <similarity treshold >" << std::endl;
    return 1;
  }

  // read the similarity threshold from argument 2
  double threshold = atof(argv[2]);
  if (threshold < 0 || threshold > 1) {
    std::cout << "Use a sensible number of iterations." << std::endl
              << "Usage: " << argv[0]
              << " <molecule-file> <similarity treshold >" << std::endl;
    return 1;
  }

  // get a molecule for each molecule entry in the provided file
  MoleculeVector mols = MoleculeFactory::getAllMolecules(std::string(argv[1]));
  unsigned int n = mols.size();

  // cluster molecules and create output
  double** M = Naomini::MakeDistanceMatrix(mols);
  std::tuple<unsigned int, unsigned int> result;
  unsigned int count = n; //Algorithmus soll theoretisch durchgefuehrt werden, bis nur noch ein Eintrag da ist. count wird mit jeder Iteration verkleinert
  std::vector<std::vector<unsigned int>> cluster;
  double greaterThan = 0; //Um nicht immer auf den gleichen kleinsten Wert zuzugreifen, wird der letzte kleinste Wert in greaterThan gespeichert
  while (count>1){
      result = Naomini::GetIndices(M, threshold, greaterThan, n);
      unsigned int low_i = std::get<0>(result);
      unsigned int low_j = std::get<1>(result);
      greaterThan = M[low_i][low_j];
      AppendToCluster(cluster, low_i, low_j);
      M=Naomini::OverwriteMatrix(M, low_i, low_j, n);
      count--;
  }
  for (unsigned int a=0; a<cluster.size();a++){
      for (unsigned int b=0; b<cluster[a].size();b++){
          std::cout << "Cluster "<< a << " containing Molecule " << cluster[a][b] <<std::endl;
    }
}
  
  
  
  return 0;
}


