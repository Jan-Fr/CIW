#pragma once

/*----------------------------------------------------------------------------*/
/**********   INCLUDES                                               **********/
/*----------------------------------------------------------------------------*/

#include <vector>
#include "Naomini/Forward.hpp"
/*----------------------------------------------------------------------------*/
/**********   NAMESPACE                                              **********/
/*----------------------------------------------------------------------------*/

namespace Naomini {

/** Calculates a non-unique, but all-encompassing set of rings
 * of the molecule. Cycles (rings) are stored as bond sets (CycleSets). Returns a
 * vector of bond sets (CyclicBondsVector).
 *
 * @brief calculates rings of a molecule
 */
  CyclicBondsVector moleculeGetRings(MoleculePtr mol);

/** Calculates a set of rings of the molecule (@see moleculeGetRings) and
 * extends every ring by the neighbor atoms that are exclusively connected to
 * the ring.
 *
 * @brief calculates rings of a molecule including one-atom ring substituents
 */
  CyclicBondsVector moleculeGetExtendedRings(MoleculePtr mol);

/** Returns rings, linker and other atoms of a molecule in their respective atom set.
 *
 * @brief calculates biconnected components
 */
  AtomSetVector moleculeGetDiverseAtomSets(MoleculePtr mol);

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
