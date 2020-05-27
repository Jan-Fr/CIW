#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <map>

#include "Naomini/Forward.hpp"
#include "Naomini/Atom.hpp"
//inlcude bonds as well now
#include "Naomini/Bond.hpp"
//FenderFroechtlingRahlf
namespace Naomini {
  //Klasse Node
  class Node {

  public:
    Node(unsigned element, AtomPtr atom1, AtomPtr atom2)
      : element(element),
        atom1(atom1),
        atom2(atom2) {}

    unsigned element;
    AtomPtr atom1;
    AtomPtr atom2;
    //gibt wieder ob elemente sich entsprechen
    bool operator==(const Node& n) const {
      return (element == n.element)
             && (atom1->getID() == n.atom1->getID())
             && (atom2->getID() == n.atom2->getID());
    }
    
    bool operator<(const Node& n) const {
      if (element != n.element) {
        return element < n.element;
      } else if (atom1->getID() != n.atom1->getID()) {
        return atom1->getID() < n.atom1->getID();
      } else {
        return atom2->getID() < n.atom2->getID();
      }
    }

    friend std::ostream& operator<<(std::ostream& os, const Node& n) {
      os << n.element << "(" << n.atom1->getID() << "," << n.atom2->getID() << ")";
      return os;
    }
  };
//Implementierung der Klasse Edge
class Edge{
  //konstruktor der klasse edge, edge hat 2 parameter:  Knoten 1 und  Knoten2; um welchen bindungstyp es sich handelt wird hier nicht miteinbezogen, da Bindungsordnungen, Stereochemie und Ladungen sowie Wasserstoffatome ignoriert werden sollen
  //access specifier ist immer public
 public: 
  Edge(Node node1, Node node2): node1(node1),node2(node2){}
  //Attribute der Klasse edge:
  Node node1;
  Node node2;
  //boolean operator to check whether node1 in fragment and structure and node 2 in fragment and structure are the same, if one of them isn't it will return false;
  //beim Vergleich einer Kante müssen die Knotenlabel verglichen werden
  bool operator==(const Edge& e) const{
    return (node1==e.node1) && (node2==e.node2);
  }
};
/** Graph Klasse zur Nutzung in BronKerbosch.hpp */

  class CompatibilityGraph {

  public:
    CompatibilityGraph(MoleculePtr mol1, MoleculePtr mol2);

    ~CompatibilityGraph();
    //funktionen des kompatibilitaetsgraphen, die in .cpp definiert werden

    size_t getNofNodes();

    std::vector <Node> getNodes() const;

    Node getNode(unsigned i) const;

    bool hasEdge(const Node& a, const Node& b) const;

  private:
/* Hier bitte Membervariablen zum Speichern von Knoten und Kanten einfuegen. */
  //zwei Vektoren, die Knoten und Kanten
  std::vector<Node> nodes;
  std::vector<Edge> edges;
  };

}
