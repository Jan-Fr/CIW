/**
 @file
 @brief
 @author    Bietz
 @date      Apr 10, 2012
 @if DoNotInclude
 Copyright ZBH  Center for Bioinformatics
 University of Hamburg
 Bundesstrasse 43, 20146 Hamburg, Germany
 ===============================================================================
 This module is part of the molecule software library Naomi,
 developed at the Center for Bioinformatics Hamburg (ZBH).
 It is not allowed to use, modify, copy or redistribute this software without
 written approval of the ZBH. For further information contact

 ZBH  Center for Bioinformatics, University of Hamburg
 Research Group for Computational Molecular Design

 Voice:  +49 (40) 42838 7350
 Fax:    +49 (40) 42838 7352
 E-Mail: info@zbh.uni-hamburg.de
 ===============================================================================
 @endif
 
 
Zeit ohne refinement:
user    23m42.642s
sys     0m0.176s

Zeit mit Refinement:
user    0m1.633s
sys     0m0.028s




 */
#include <iomanip>
#include <stdio.h>

#include "Ullmann.hpp"
#include "Naomini/Forward.hpp"
#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Helpers.hpp"
#include "Naomini/MoleculeDrawer.hpp"

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE NAOMINI                                      **********/
/*----------------------------------------------------------------------------*/
namespace Naomini {

/*----------------------------------------------------------------------------*/
/**********   NAMELESS NAMESPACE                                     **********/
/*----------------------------------------------------------------------------*/
  namespace {

/*------------*/
/** HELPERS  **/
/*------------*/

    bool atoms_are_adjacent(AtomPtr a, AtomPtr b) {
      for (AtomPtr n : a->getNeighborAtoms()) {
        if (b == n) {
          return true;
        }
      }
      return false;
    }

/*----------------------------------------------------------------------------*/

    BoolMatrix get_adjacency_matrix(const AtomVector& atoms) {
      unsigned nofAtoms = atoms.size();
      BoolMatrix adj(nofAtoms, nofAtoms);
      for (unsigned i = 0; i < nofAtoms - 1; ++i) {
        adj(i, i) = false;
        for (unsigned j = i + 1; j < nofAtoms; ++j) {
          adj(i, j) = adj(j, i) = atoms_are_adjacent(atoms.at(i), atoms.at(j));
        }
      }
      return adj;
    }

/*----------------------------------------------------------------------------*/

    unsigned nofHeavyNeighbors(AtomPtr atom) {
      unsigned count = 0;
      for (AtomPtr n : atom->getNeighborAtoms()) {
        if (!atomIsHydrogen(n)) count++;
      }
      return count;
    }

/*----------------------------------------------------------------------------*/
/**********   END NAMELESS NAMESPACE                                 **********/
/*----------------------------------------------------------------------------*/
  }


/*---------------*/
/** CONSTRUCTOR **/
/*---------------*/

  Ullmann::Ullmann(MoleculePtr substructure, MoleculePtr molecule)
    : m_Molecule(molecule) {
    for (AtomPtr sAtom : substructure->getAtoms()) {
      if (!atomIsHydrogen(sAtom)) {
        m_SubstructureAtoms.push_back(sAtom);
      }
    }
    for (AtomPtr mAtom : molecule->getAtoms()) {
      if (!atomIsHydrogen(mAtom)) {
        m_MoleculeAtoms.push_back(mAtom);
      }
    }

    m_MoleculeAdjacencyMatrix = get_adjacency_matrix(m_MoleculeAtoms);
    m_SubstructureAdjacencyMatrix = get_adjacency_matrix(m_SubstructureAtoms);
  }

/*------------*/
/** METHODS  **/
/*------------*/

  void Ullmann::run() {
    BoolMatrix initMatrix = getInitialMatrix();
    printMatrix(initMatrix, 0, "Initialmatrix");

    BoolVector F(m_MoleculeAtoms.size(), false);
    enumerateSubgraph(initMatrix, 0, F);
  }

/*----------------------------------------------------------------------------*/

  bool Ullmann::checkSubgraph(const BoolMatrix& M) {
    using namespace boost::numeric::ublas;

    BoolMatrix matchingAdj = prod(M, trans(prod(M, m_MoleculeAdjacencyMatrix)));

    for (unsigned i = 0; i < m_SubstructureAdjacencyMatrix.size1() - 1; ++i) {
      for (unsigned j = i + 1; j < m_SubstructureAdjacencyMatrix.size2(); ++j) {
        if (m_SubstructureAdjacencyMatrix(i, j) && !matchingAdj(i, j)) {
          return false;
        }
      }
    }
    return true;
  }

/*----------------------------------------------------------------------------*/

  void Ullmann::drawerSubstructureMatching(const BoolMatrix& M) {
    AtomVector substructure;

    for (unsigned i = 0; i < M.size1(); ++i) {
      for (unsigned j = 0; j < M.size2(); ++j) {
        if (M(i, j)) {
          substructure.push_back(m_MoleculeAtoms.at(j));
          break;
        }
      }
    }

    MoleculeDrawer drawer(m_Molecule);
    drawer.markSubstructure(substructure, MoleculeDrawer::GREEN);

  }

/*----------------------------------------------------------------------------*/

  void Ullmann::enumerateSubgraph(BoolMatrix& M, unsigned k, BoolVector F) {

    if (k == m_SubstructureAtoms.size()) {
      if (checkSubgraph(M)) {
        drawerSubstructureMatching(M);
      }
    } else {
      for (unsigned l = 0; l < M.size2(); ++l) {
        if (M(k, l) && !F.at(l)) {
          BoolMatrix copyOfM = M;
          for (unsigned j = 0; j < M.size2(); ++j) {
            M(k, j) = false;
          }
          M(k, l) = true;
          F.at(l) = true;
//        printMatrix(M, k+1, "");
          if(Ullmann::refine(M,k)){
            enumerateSubgraph(M, k + 1, F);
          }
          F.at(l) = false;
          M = copyOfM;
        }
      }
    }
  }
  
  
  bool Ullmann::refine(BoolMatrix& M, unsigned k){
    bool change;
    //Ullmann Kriterium
    bool crit;
    bool neighbor;
    bool possibleMatch;
    do{
    change = false;
    // i iteriert durch Zeilen, j durch Spalten von M
    for (unsigned i= k+1; i<M.size1(); ++i){
      for (unsigned j=0; j<M.size2(); ++j){
        if(M(i,j)){
          crit = true;
          // Nachbaratome der Substruktur und des Molekuls in AtomVektoren speichern
          AtomVector subH = m_SubstructureAtoms.at(i)->getNeighborAtoms();
          AtomVector sub;
          AtomVector molH = m_MoleculeAtoms.at(j)->getNeighborAtoms();
          AtomVector mol;
          // Nur nicht-H-Nachbaratome sollen betrachtet werden
          for (unsigned s=0; s<subH.size(); s++){
            if (!atomIsHydrogen(subH.at(s))){
              sub.push_back(subH.at(s));    
            }
          }
          // Nur nicht-H-Nachbaratome sollen betrachtet werden
          for (unsigned t=0; t<molH.size(); t++){
            if (!atomIsHydrogen(molH.at(t))){
              mol.push_back(molH.at(t));    
            }
          }
          // Iteration durch die Nachbaratome des aktuellen Atoms der Substruktur und des Molekuels
          for (unsigned l = 0; l<sub.size();++l){
            neighbor = false;
            for (unsigned m = 0; m<mol.size();m++){
                //getID gibt den Index eines Atoms zurueck
                //Ueberpruefung die Nachbaratome der Substruktur mit einem der Nachbaratome des Molekuels matcht
              if (M(sub.at(l)->getID(),mol.at(m)->getID())){
                neighbor = true;
                break;
              }
            }
            // Falls es keine passenden Nachbarn gibt, ist das Ullmannkriterium nicht erfuellt und es wird abgebrochen
            if(!neighbor){
                crit = false;
                break;
            }
          }
          //Wenn Ullmannkriterium nicht erfuellt ist, wird die bool Zuweisung in M auf 0 gesetzt
          //Permutationsaenderung(change) wird auf true gesetzt
          if(!crit){
            M(i,j) = 0;
            change = true;
            possibleMatch = false;
            for (unsigned h=0; h<M.size2();h++){
              if(M(i,h)){
                possibleMatch = true;
              }
            }
            if(!possibleMatch){
              return false;
            }
          }
        }
      }
    }
  }while(change);
  return true;
      
  }

/*----------------------------------------------------------------------------*/

  BoolMatrix Ullmann::getInitialMatrix() {
    unsigned nofSubstrAtoms = m_SubstructureAtoms.size();
    unsigned nofMolAtoms = m_MoleculeAtoms.size();

    BoolMatrix initMatrix(nofSubstrAtoms, nofMolAtoms);
    for (unsigned i = 0; i < nofSubstrAtoms; ++i) {
      AtomPtr subAtom = m_SubstructureAtoms.at(i);
      for (unsigned j = 0; j < nofMolAtoms; ++j) {
        AtomPtr molAtom = m_MoleculeAtoms.at(j);
        if ((subAtom->getAtomicNumber() == molAtom->getAtomicNumber()) &&
            nofHeavyNeighbors(subAtom) <= nofHeavyNeighbors(molAtom)) {
          initMatrix(i, j) = true;
        } else {
          initMatrix(i, j) = false;
        }
      }
    }
    return initMatrix;
  }

/*----------------------------------------------------------------------------*/

  void Ullmann::printMatrix(
    BoolMatrix& M,
    unsigned k,
    const std::string& name)
  {
    std::cout << name << std::endl;
    std::cout << " " << std::setw(3) << k << " ║";
    for (unsigned j = 0; j < M.size2(); j++) {
      std::cout << std::setw(4) << std::setprecision(2) << m_MoleculeAtoms.at(j)->getName() << " │";
    }
    std::cout << "\n═════╬═════";
    for (unsigned j = 1; j < M.size2(); j++) {
      std::cout << "╪═════";
    }
    std::cout << "╡\n";
    for (unsigned i = 0; i < M.size1(); i++) {
      std::cout << " " << std::setw(3) << m_SubstructureAtoms.at(i)->getName() << " ║";
      for (unsigned j = 0; j < M.size2(); j++) {
        if (M(i, j)) {
          printf("  1  ");
        } else {
          printf("     ");
        }
        std::cout << "│";
      }
      std::cout << "\n─────╫─────";
      for (unsigned j = 1; j < M.size2(); j++) {
        std::cout << "┼─────";
      }
      std::cout << "┤\n";
    }
    std::cout << std::endl;
  }


/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
