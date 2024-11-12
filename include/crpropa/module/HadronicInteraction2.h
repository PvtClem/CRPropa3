#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include <fstream>
#include <cmath>

#include "crpropa/Module.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class HadronicInteraction
 @brief Hadronic interactions for some cosmic rays


*/
class HadronicInteraction: public Module {
private:

  double Density;
	bool haveElectrons;
  bool havePhotons;
  bool haveNeutrinos;
	double limit;
	double thinning;
	std::string interactionTag = "HI";

	// tabulated interaction rate 1/lambda(E)
  std::vector<double> tabEnergy;  //!< hadron energy in [GeV]
	std::vector<double> tabCS;  //!< cross-section

	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< hadron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate

public:
	/** Constructor
   @param Density
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param havePhotons	if true, add secondary photons as candidates
	 @param haveNeutrinos	if true, add secondary neutrinos as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 */
	HadronicInteraction(double Density = 0.0, bool haveElectrons = false, bool havePhotons = false, bool haveNeutrinos = false, double thinning = 0, double limit = 0.1);

  /**
  @param Density
  */
  void setGasDensity(double Density);

	// decide if secondary electrons are added to the simulation
	void setHaveElectrons(bool haveElectrons);

  	// decide if secondary photons are added to the simulation
	void setHavePhotons(bool havePhotons);

  	// decide if secondary neutrinos are added to the simulation
	void setHaveNeutrinos(bool haveNeutrinos);

	/** limit the step to a fraction of the mean free path
	 @param limit	fraction of the mean free path, should be between 0 and 1
	*/
	void setLimit(double limit);

	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);

	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;

	void getCS(std::string filename) ;

	void process(Candidate *candidate)  ;
	void performInteraction(Candidate *candidate) const ;

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H
