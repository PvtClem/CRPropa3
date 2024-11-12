#ifndef CRPROPA_HADRONIC_H
#define CRPROPA_HADRONIC_H

#include <fstream>
#include <cmath>

#include "crpropa/Module.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class Hadronic
 @brief Hadronic interactions for some cosmic rays


*/
class Hadronic: public Module {
private:



public:
	/** Constructor

	 */
	Hadronic();


	void process(Candidate *candidate) const ;

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_HADRONIC_H
