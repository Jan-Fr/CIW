#include "RingsAndBiconnectedComponents.hpp"

#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Forward.hpp"
#include "Naomini/Helpers.hpp"
#include "Naomini/Molecule.hpp"

using namespace Naomini;

#include <vector>
#include <algorithm>
#include <map>
#include "Naomini/MoleculeFactory.hpp"
#include "Naomini/MoleculeDrawer.hpp"


//Implementierung der Tiefensuche
void dfs(std::map<int,bool> &visited,AtomPtr currentAtom,int index, BondPtr e);
void dfs(std::map<int,bool> &visited,AtomPtr currentAtom,int index, BondPtr e){
  //aktueller Knoten wird als besucht markiert
  visited[index] = true;
  //Iteration durch Nachbaratome
  for (AtomPtr a: currentAtom->getNeighborAtoms()){
    //Ein pair aus currentAtom und dem betrachteten Nachbaratom a bilden
    //bzw zwei pairs, da (a, current) und (current, a) moeglich sind
    std::pair<AtomPtr, AtomPtr> connectedAtoms1;
    std::pair<AtomPtr, AtomPtr> connectedAtoms2;
    connectedAtoms1.first = a;
    connectedAtoms1.second = currentAtom;
    connectedAtoms2.first = currentAtom;
    connectedAtoms2.second = a;
    //Die Atome der entfernten Kante e werden mit den pairs verglichen
    //Bindungen der entfernten Kante sollen ignoriert werden
    if (e->getAtoms() != connectedAtoms1 && e->getAtoms() != connectedAtoms2){
      //Verschieben des Index auf den Index von a
      index=a->getID();
      //dfs soll nur bei noch nicht besuchten Atomen erneut aufgerufen werden
      if(!visited[index]){
        dfs(visited, a, index, e);
      }
    }
  }
}

CyclicBondsVector Naomini::moleculeGetRings(MoleculePtr mol) {
  // Insert your code here but do not use this function for your solution!
  //CyclicBondsVector rings = getCyclicSetsOfMolecule(mol);
  
  AtomVector atoms = mol->getAtoms();
  CyclicBondsVector rings;
  //Alle Kanten betrachten, eine zyklische Bindung ist gegeben, wenn trotz entfernen der Kante, alle Atome durch dfs erreicht werden koennen
  for(BondPtr e: mol->getBonds()){
    //visited map mit false fuer jedes Atom initialisieren
    std::map<int,bool> visited;
    for (AtomPtr a: atoms){
      visited[a->getID()] = false;
    }
    AtomPtr currentAtom = atoms.at(0);
    int index = 0;
    dfs(visited,currentAtom, index, e);
    //nach dfs wird visited geprueft. Sollte einmal false vorkommen wird isCyclic auf false gesetzt
    bool isCyclic = true;
    for (AtomPtr a: atoms){
      if(visited[a->getID()]== false){
        std::cout << a->getID() <<std::endl;
        isCyclic = false;
      }
    }
    
    //Wenn alle Elemente in visited true sind, ist auch isCyclic true
    if(isCyclic==true){
      //Keine Ahnung ob es noch einen weniger umstaendlichen weg gibt ein Element in rings einzufuegen, aber das hat geklappt
      std::set<BondPtr> add;
      add.insert(e);
      rings.push_back(add);
    }
  }
  return rings;
}


/*----------------------------------------------------------------------------*/

CyclicBondsVector Naomini::moleculeGetExtendedRings(MoleculePtr mol) {
  CyclicBondsVector rings;
  rings = moleculeGetRings(mol);
  //rings ist ein vector von sets
  //iteration durch vektor
  for (unsigned int i=0; i<rings.size(); ++i){
      //iteration durch set
    for(std::set<BondPtr> :: iterator it = rings[i].begin(); it != rings[i].end(); ++it){
      //Die Bindung liegt in *it
      BondPtr bond = *it;
      //Beteiligte Atome der Bindung extrahieren und in Vektor speichern
      AtomPair ap = bond->getAtoms();
      AtomVector av;
      av.push_back(ap.first);
      av.push_back(ap.second);
      //Iteration durch die zwei an der Bindung beteiligten Atome
      for(AtomPtr a: av){
        //Nachbarn der Atome betrachten und ueber sie iterieren
        AtomVector neighbors = a->getNeighborAtoms();
        for (AtomPtr n: neighbors){
          //Bedingungen: Endstaendig und HeavyAtom
          if (n->getNofBonds() == 1 && n->getElementSymbol()!="H"){
            //Einfuegen der Bindung in rings (Auch hier geht das bestimmt noch schoener)
            std::set<BondPtr> add;
            for(BondPtr bondToAdd: n->getBonds()){
              add.insert(bondToAdd);
              rings.push_back(add);
            }
          }
        }
      }
    }
  }
  return rings;
}

/*----------------------------------------------------------------------------*/

AtomSetVector Naomini::moleculeGetDiverseAtomSets(MoleculePtr mol) {
  AtomSetVector atomSets;
  throw "Insert your code in RingsAndBiconnectedComponents.cpp";
  return atomSets;
}

