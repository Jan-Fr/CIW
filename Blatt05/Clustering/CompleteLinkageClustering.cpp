#include "CompleteLinkageClustering.hpp"
#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeDrawer.hpp"
#include "Naomini/MoleculeFactory.hpp"



using namespace Naomini;

double** Naomini::MakeDistanceMatrix(MoleculeVector mols){
    int n = mols.size();
    double** M;
    M = new double*[n];
    for (unsigned i=0; i<mols.size();i++){
        M[i] = new double[n];
        for (unsigned j=0; j<mols.size();j++){
            ECFP a = getECFPOfMolecule(mols.at(i), 2);
            
            ECFP b = getECFPOfMolecule(mols.at(j), 2);
            
            double tani = getTanimotoCoefficientOfECFPs(a,b);
            
            M[i][j] = tani;
        }
    }
    return M;
}


std::tuple<unsigned int, unsigned int> Naomini::GetIndices(double** M, double sim, double greaterThan, unsigned int n){
    unsigned int low_i = 0;
    unsigned int low_j = 0;
    double min = M[0][0];
    for (unsigned int i=0; i<n; i++){
        for (unsigned int j=0; j<n; j++){
            if(M[i][j]<min && M[i][j]>greaterThan){
                min = M[i][j];
                low_i=i;
                low_j=j;
            }
        }
    }
    //Tupel der Indices des kleinsten Eintrags erstellen
    //Wenn threshhold ueberschritten wurde, wird NULL zurueckgegeben
    std::tuple<unsigned int, unsigned int> result = std::make_tuple(NULL,NULL);
    if(M[low_i][low_j]<=sim){
        return std::make_tuple(low_i,low_j);
    }
    return result;
}



double** Naomini::OverwriteMatrix(double** M, unsigned int low_i, unsigned int low_j, unsigned int n){
    /*geht die komplette Matrix durch. Wenn in Zeile/Spalte der Molekuele mit
    kleinstem Abstand: Maxima der Abstaende zu den anderen Molekuelen verwenden */
    for (unsigned int i=0; i<n; i++){
        for (unsigned int j=0; j<n;j++){
            if (i==low_i || i==low_j){
                M[i][j] = std::max(M[low_i][j],M[low_j][j]);
            }
            else if (j==low_j || j == low_i){
                M[i][j] = std::max(M[i][low_i],M[i][low_j]);
            }
        }
    }
    return M;
}

void Naomini::AppendToCluster(std::vector<std::vector<unsigned int>> &cluster, unsigned int low_i, unsigned int low_j){
    bool appended = false;
    unsigned int n = cluster.size();
    // cluster ist in der ersten Iteration noch leer
    if (cluster.empty() == false){
        for (unsigned int i=0; i<n; i++){
            for (unsigned int j=0; j<cluster[i].size(); j++){
                /*In cluster sind die Indices der kleinsten Distanzen gespeichert.
                Wenn ein neues Indices-Paar ein Cluster bilden soll, wird hier ueberprueft, ob es einer der Indices schon in einem Cluster vorhanden ist. Falls ja, wird der andere Index dem gleichen Cluster hinyugefuegt.*/
                if (cluster[i][j] == low_i){
                    cluster[i].push_back(low_j);
                    appended = true;
                    break;
                }
                else if (cluster[i][j] == low_j){
                    cluster[i].push_back(low_i);
                    appended = true;
                    break;
                }
            }
            if (appended == true){
                break;
            }
        }
        //Falls beide Indices noch keinem Cluster angehoeren, bilden sie ein neues Cluster
        if (appended == false){
            cluster.push_back({low_i, low_j});
        }
    }
    else{
        //Erste Iteration
        cluster = {{low_i,low_j}};
    }
}
