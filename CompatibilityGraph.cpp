#include "CompatibilityGraph.hpp"

#include <stdexcept>

#include "Naomini/Forward.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Bond.hpp"
//FenderFroechtlingRahlf
namespace Naomini {

/* Helper Function ************************************************************/
  namespace {
    //boolean moleculeHasBond -->iteriert durch die bonds, speichert alle Atompaare in atoms, wenn atompaare existieren existiert ein bond und true wird zurückgegeben
    bool moleculeHasBond(MoleculePtr mol, AtomPtr atom1, AtomPtr atom2) {
      for (BondPtr bond : mol->getBonds()) {
        AtomPair atoms = bond->getAtoms();
        if ((atoms.first == atom1 && atoms.second == atom2) ||
            (atoms.first == atom2 && atoms.second == atom1)) {
          return true;
        }
      }
      return false;
    }

  } // end anonymouse namespace
/* Helper Function END ********************************************************/

/* Compatibility Class Functions **********************************************/

  CompatibilityGraph::CompatibilityGraph(MoleculePtr mol1, MoleculePtr mol2) {
/* Hier bitte den Code zur Initialisierung des Graphen einfuegen. */
    //throw std::runtime_error("Constructor not implemented");
    //Iterieren durch Knoten des Moleküls 1 
    for(AtomPtr atom1: mol1->getAtoms()){
      //Ignorieren von Wasserstoffen, dann einfach continue zum nächsten Element
      if (atom1->getElementName()=="Hydrogen"){
        continue;
      }
      //Iterieren durch Knoten des Moleküls 2
      for(AtomPtr atom2:mol2->getAtoms()){
        //wieder wenn Wasserstoff-->continue
        if(atom2->getElementName()=="Hydrogen"){
          continue;
        }  
        //wenn atom1 und atom2 das gleiche element haben, dann werden sie zu dem vector nodes hinzugefügt
        if(atom1->getElementName()==atom2->getElementName()){
          //Zufügen zu nodes: ein Node besteht aus der ElementID (hier für beide Atome gleich) und zwei AtomPointern
          nodes.push_back(Node(atom1->getAtomicNumber(), atom1, atom2));
        }
      }
    }
    //Iterieren durch die zuvor erstellten nodes, die das gleiche Element hatten
    for (Node node1:nodes){
      for (Node node2:nodes){
        //wir wollen nicht den gleichen Knoten als node1 und node2 haben
        if(node1==node2){
          continue;
        }
        //wenn zwischen der knoten in mol1 und zwischen den knoten in mol2 eine kante besteht, speichern der Kante im Vector edges
        if(moleculeHasBond(mol1,node1.atom1,node2.atom1)==moleculeHasBond(mol2,node1.atom2,node2.atom2)){
          //Edge besteht aus Knoten 1 und Knoten 2
          edges.push_back(Edge(node1,node2));
        }
      }
    }
    
  }

  CompatibilityGraph::~CompatibilityGraph() {}

  size_t CompatibilityGraph::getNofNodes() {
/* Anzahl der Knoten */
    //gibt grösse des vectors nodes wieder
    return nodes.size();
    //throw std::runtime_error("getNofNodes not implemented");
  }

  std::vector <Node> CompatibilityGraph::getNodes() const {
/* Hier bitte den Code zur Abfrage der Knoten eingeben. */
    //gibt vector nodes wieder
    return nodes;
    //throw std::runtime_error("getNodes not implemented");
  }

  Node CompatibilityGraph::getNode(unsigned i) const {
/* Hier bitte den Code zur Abfrage eines Knotens eingeben. */
    //-->gib Knoten an Stelle i in nodes wieder
    return nodes[i];
    //throw std::runtime_error("getNode not implemented");
  }

  bool CompatibilityGraph::hasEdge(const Node& a, const Node& b) const {
/* Hier bitte den Code zur Abfrage der Existenz einer Kante einfuegen. */
    return (find(edges.begin(),edges.end(),Edge(a,b))!=edges.end());
    //throw std::runtime_error("hasEdge not implemented");
  }

}

