#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/massDistribution/ConstantDensity.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

HadronicInteraction::HadronicInteraction( double Density, bool haveElectrons, bool havePhotons, bool haveNeutrinos, double thinning, double limit) {
  setGasDensity(Density);
	setHaveElectrons(haveElectrons);
	setHavePhotons(havePhotons);
	setHaveNeutrinos(haveNeutrinos);
	setLimit(limit);
	setThinning(thinning);
}


void HadronicInteraction::setGasDensity(double Density) {
	this->Density = Density;
}

void HadronicInteraction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void HadronicInteraction::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void HadronicInteraction::setHaveNeutrinos(bool haveNeutrinos) {
	this->haveNeutrinos = haveNeutrinos;
}

void HadronicInteraction::setLimit(double limit) {
	this->limit = limit;
}

void HadronicInteraction::setThinning(double thinning) {
	this->thinning = thinning;
}

void HadronicInteraction::getCS(std::string filename) {

  std::ifstream infile(filename.c_str());
  tabCS.clear(); // clear the tables
  tabEnergy.clear();


	if (!infile.good())
		throw std::runtime_error("HadronicInteraction: could not open file " + filename);

	// clear previously loaded interaction rates


	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(a * GeV);
				tabCS.push_back(b * barn);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n'); //to ignore the rest of the line
	}
	infile.close();


}


void HadronicInteraction::performInteraction(Candidate *candidate) const {

  int id = candidate->current.getId();

  std::string filename;
  filename = "./../../../../Tables_HI/cum_el_p_p04.txt";
  std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("HadronicInteraction: could not open file " + filename);

	// clear previously loaded interaction rates, and then read our 2d data table


  std::vector<std::vector<double>> tabdCS; // differential cross-section
  std::vector<double> tabE_sec;  //!< hadron energy in [GeV]


  bool firstLine = true;
  double value;
  double firstValue;
  double E = candidate->current.getEnergy();

	    while (infile.good()) {
        if (infile.peek() != '#') {  // Ignore lines starting with '#'
            std::string line;
            if (std::getline(infile, line)) {
                std::istringstream stream(line);
                // stream >> firstValue; // we don't want the first value

                if (firstLine) {
                    while (stream >> value){
                    tabE_sec.push_back(value * GeV);  // Convert and store each remaining value
                }

                    firstLine = false;
                }

                // Store the rest of the line into a row vector
                std::vector<double> row;
                while (stream >> value) {
                    row.push_back(value * barn);  // Convert and store each remaining value
                }

                // Add row to tabdCS
                if (!row.empty()) {   // if row is not empty
                    tabdCS.push_back(row);
                }
            }
        }
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Skip to the next line
    }

    infile.close();

  std::vector<double> dCS_interp;

  for (size_t i = 0; i < tabdCS[0].size(); ++i) {
        dCS_interp.push_back(interpolate(E, tabdCS[i], tabE_sec)) ;
    }


	Random &random = Random::instance();

  double random_number = random.rand(); // between 0 and 1

	double E_sec = interpolate(random_number, dCS_interp, tabE_sec);

	double f = E_sec / E;

	if (haveElectrons) {
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(11, E_sec, pos, w, interactionTag);
		}

	}
	// Update the primary particle energy.
	// This is done after adding the secondaries to correctly set the secondaries parent
	candidate->current.setEnergy(E - E_sec );

  tabE_sec.clear();
	tabdCS.clear();
}

void HadronicInteraction::process(Candidate *candidate)  {
  // check if electron / positron


  int id = candidate->current.getId();
if ((abs(id) != 1000010010) && (abs(id) != 1000020040) && (abs(id) != 1000060120) && (abs(id) != 1000130260) && (abs(id) != 1000260520 )) // we don't want particles that are not in the tables (H, He, C, Al, Fe)
		return;

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();

if (abs(id) == 1000010010){
  getCS("./../../../../Tables_HI/el_p_p04.txt");
}

	// check if in tabulated energy range : suppose we have already read the data files ; we are comparing the energy of the primary here
	if ((E < tabEnergy.front()) or (E > tabEnergy.back()))
		return;

  double cross_section = interpolate(E, tabEnergy, tabCS);
  Vector3d position = candidate->current.getPosition();
  double inverse_mfp = cross_section * Density / barn ;

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / inverse_mfp;  //The negative logarithm of a uniform random variable transforms it into an exponentially distributed random variable
		// check for interaction; if it doesn't occur, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / inverse_mfp);
			return;
		}
		performInteraction(candidate);
		step -= randDistance;
	} while (step > 0.);
}


void HadronicInteraction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string HadronicInteraction::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
