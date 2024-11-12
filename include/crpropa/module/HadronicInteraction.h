#ifndef CRPROPA_HadronicInteraction_H
#define CRPROPA_HadronicInteraction_H

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


  double H_density;
  double He_density;
	bool haveElectrons;
  bool havePhotons;
  bool haveNeutrinos;
	double limit;
	double thinning;
	std::string interactionTag = "HadronicInteraction";

  std::map<int, std::string> part_numbering;
  std::vector<std::string> Filenames;
  std::map<std::string, std::vector<std::vector<double>>> dict; // E primary and cross section, at a gievn E primary (ie integrated over Esec)
  std::map<std::string, std::vector<double>> dictE_sec; // list of E sec
  std::map<std::string, std::vector<double>> dictE_prim; // list of E prim that are in cumulative data
  std::map<std::string, std::vector<std::vector<double>>> dict_cumDCS; // cumulative differential cross-section


public:
	/** Constructor
   @param H_density; H density in SI
   @param He_density; He density in SI
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param havePhotons	if true, add secondary photons as candidates
	 @param haveNeutrinos	if true, add secondary neutrinos as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 */
	HadronicInteraction(double H_density = 0.0, double He_density = 0.0, bool haveElectrons = false, bool havePhotons = false, bool haveNeutrinos = false, double thinning = 0, double limit = 0.1);

  /**
  @param Density
  */
  void setGasDensity(double H_density, double He_density);

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

  void setTables();

  std::vector<double> setInteraction(Candidate *candidate) const;

	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;


	void process(Candidate *candidate) const ;
	void performInteraction(Candidate *candidate, double secondary) const ;

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_HadronicInteraction_H

