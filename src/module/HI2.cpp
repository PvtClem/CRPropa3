#include "crpropa/module/HI2.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/massDistribution/ConstantDensity.h"

#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>


#include <iostream>



namespace crpropa {

static const double mec2 = mass_electron * c_squared * eV / 1e-9; // in GeV
static const double mpc2 = mass_proton * c_squared * eV / 1e-9;


HI::HI( double density, bool haveElectrons, bool havePhotons, bool haveNeutrinos, double thinning, double limit) {
  setGasDensity(density); // this one will probably disappear in order to be set from the simulation
	setHaveElectrons(haveElectrons);
	setHavePhotons(havePhotons);
	setHaveNeutrinos(haveNeutrinos);
	setLimit(limit);
	setThinning(thinning); // functions run to initialize the HI class
}


void HI::setGasDensity(double density) {
	this-> density = density;
}

void HI::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void HI::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void HI::setHaveNeutrinos(bool haveNeutrinos) {
	this->haveNeutrinos = haveNeutrinos;
}

void HI::setLimit(double limit) {
	this->limit = limit;
}

void HI::setThinning(double thinning) {
	this->thinning = thinning;
}

std::vector<double> HI::setInteraction(Candidate *candidate) const { // given a particle, get the primary total energies back, and use it (and the cross-sections) to decide what secondary will be produced (provided an interaction occurs). Then we need to decide whether the interaction occurs or not (given the total cross-section computed here, in process), and then perform the interaction (in performInteraction)

  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) / eV * 1e-9; // convert energy from J to GeV

  double CS_tot = 0;
  double secondary = 0;

  std::vector<double> tabCS_E; // table of cross sections (for different processes) at a given energy : used to decide which process we perform
  std::string name;

  std::vector<double> tabEnergy;
  std::vector<std::vector<double>> tabCS(5); // differential cross-section, to have at least a sisze of 5 ?

  std::vector<std::string> tab_filenames;

  //std::cout << "Hello, World!" << std::endl;


  /** to get the current working directory
  char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        perror("getcwd() error");
    }
  **/

  if (abs(id) == 1000010010){ // if we have protons
     tab_filenames = {"el_p_p04", "pos_p_p04", "prot_p_p04", "aprot_p_p04", "gam_p_p04"}; // suppose for now that the gas is only made of HI
  }

  else if (abs(id) == 1000020040){ // if we have protons
     tab_filenames = {"el_He_p04", "pos_He_p04", "prot_He_p04", "aprot_He_p04", "gam_He_p04"};
  }

  else if (abs(id) == 1000060120){ // if we have protons
     tab_filenames = {"el_C_p04", "pos_C_p04", "prot_C_p04", "aprot_C_p04", "gam_C_p04"};
  }

  else if (abs(id) == 1000130260){ // if we have protons
     tab_filenames = {"el_Al_p04", "pos_Al_p04", "prot_Al_p04", "aprot_Al_p04", "gam_Al_p04"};
  }

  else if (abs(id) == 1000260520){ // if we have protons
    tab_filenames = {"el_Fe_p04", "pos_Fe_p04", "prot_Fe_p04", "aprot_Fe_p04", "gam_Fe_p04"};
  }


  tab_filenames = {"el_p_p04", "pos_p_p04", "prot_p_p04", "aprot_p_p04", "gam_p_p04"}; // we only have protons for now, both on the CR side and the

  for (int i = 0; i < 5; i++){

    name = "./Tables_HI/" + tab_filenames[i] + ".txt";

    std::ifstream infile(name.c_str());

    //std::cout << "i: ";
    //std::cout << i << std::endl;


	  if (!infile.good()){
      std::cout << name << std::endl;
		  throw std::runtime_error("HI: could not open file " + name);
      //std::cout << "not good" << std::endl;
      }


    if (i==0){
      //std::cout << "there" << std::endl;
      while (infile.good()) {
		  if (infile.peek() != '#') {
        //std::cout << "in the loop" << std::endl;
			  double a, b;
			  infile >> a >> b;
			  if (infile) {
          //std::cout << "in the if" << std::endl;
				  tabEnergy.push_back(a );
				  tabCS[i].push_back(b );
          //std::cout << "passed it" << std::endl;
			  }
		  }
		  infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n'); //to ignore the rest of the line
	  }
    }
    else{ // we already have the primaries energy stored, don't need them twice'
      //std::cout << "2!" << std::endl;
      //std::cout << i << std::endl;
	    while (infile.good()) {
		    if (infile.peek() != '#') {
			    double a, b;
			    infile >> a >> b;
			    if (infile) {
				    tabCS[i].push_back(b);
			    }
		    }
		    infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n'); //to ignore the rest of the line
	    }
    }
  infile.close();
  tabCS_E.push_back(interpolate(E, tabEnergy, tabCS[i])); // get the CS at the energy of the primary interacting particle ; interpolate(x, tabX, tabY)
  CS_tot += tabCS_E[i]; // -1 is not allowed in C++
  }

  Random &random = Random::instance();
  double random_number = random.rand() * CS_tot;
  double summed_CS = tabCS_E[0];

  while (summed_CS <= random_number){
    secondary += 1;
    summed_CS += tabCS_E[secondary];
  }

  //std::cout << CS_tot ;
  //std::cout << "  " ;
  //std::cout << secondary;
  return std::vector<double> {CS_tot, secondary};
}


void HI::performInteraction(Candidate *candidate, double secondary) const {

  //std::cout << "performInteraction" << std::endl;

  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) / eV * 1e-9;

  std::vector<std::string> tab_filenames;

  if (abs(id) == 1000010010){ // if we have protons
     tab_filenames = {"el_p_p04", "pos_p_p04", "prot_p_p04", "aprot_p_p04", "gam_p_p04"}; // suppose for now that the gas is only made of HI
  }

  else if (abs(id) == 1000020040){ // if we have protons
     tab_filenames = {"el_He_p04", "pos_He_p04", "prot_He_p04", "aprot_He_p04", "gam_He_p04"};
  }

  else if (abs(id) == 1000060120){ // if we have protons
     tab_filenames = {"el_C_p04", "pos_C_p04", "prot_C_p04", "aprot_C_p04", "gam_C_p04"};
  }

  else if (abs(id) == 1000130260){ // if we have protons
     tab_filenames = {"el_Al_p04", "pos_Al_p04", "prot_Al_p04", "aprot_Al_p04", "gam_Al_p04"};
  }

  else if (abs(id) == 1000260520){ // if we have protons
    tab_filenames = {"el_Fe_p04", "pos_Fe_p04", "prot_Fe_p04", "aprot_Fe_p04", "gam_Fe_p04"};
  }

  tab_filenames = {"el_p_p04", "pos_p_p04", "prot_p_p04", "aprot_p_p04", "gam_p_p04"};
  std::string filename;
  filename = "./Tables_HI/cum_" + tab_filenames[secondary] + ".txt";


  std::ifstream infile(filename.c_str());

  //std::cout << "File read" << std::endl;

	if (!infile.good())
		throw std::runtime_error("HI: could not open file " + filename);

	// clear previously loaded interaction rates, and then read our 2d data table


  std::vector<std::vector <double>> tabdCS; // differential cross-section
  std::vector<double> tabE_sec;  //!< hadron kinetic energy in [GeV]
  std::vector<double> tabE_primary;


  bool firstLine = true;
  double value;
  double firstValue;

	    while (infile.good()) {
         //std::cout << "entering while" << std::endl;
        if (infile.peek() != '#') {  // Ignore lines starting with '#'
            std::string line;
            if (std::getline(infile, line)) {
                std::istringstream stream(line);
                stream >> firstValue;

                if (firstLine) {
                    while (stream >> value){
                    tabE_sec.push_back(value );  //  store each remaining value
                     //std::cout << "end of second while" << std::endl;
                }

                    firstLine = false;
                }

                else{
                  tabE_primary.push_back(firstValue); // don't want the first 0, corresponding to the kinetic energy line

                  // Store the rest of the line into a row vector
                  std::vector<double> row;
                  while (stream >> value) {
                    row.push_back(value );  // store each remaining value ! Careful, it's mbarn
                    //std::cout << "end of third while" << std::endl;
                  }

                  // Add row to tabdCS
                  if (!row.empty()) {   // if row is not empty
                    tabdCS.push_back(row);
                    //std::cout << "added to tabdCS" << std::endl;
                }
                }
            }
        }
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Skip to the next line
    }

    infile.close();


  std::vector<double> dCS_interp;

  for (size_t i = 0; i < tabdCS[0].size(); ++i) {
        std::vector<double> result;
        std::transform(tabdCS.begin(), tabdCS.end(), std::back_inserter(result),[i](const std::vector<double>& subvec) {return subvec[i];}); // extract ith element of each subvector of tabdCS
        dCS_interp.push_back(interpolate(E, tabE_primary, result)) ; // would need to interpolate en tabE here, and we interpolate with something like tabdCS(:)(i), where i is the column corresponding tot a given Esec
    }
  //std::cout << "passed the for" << std::endl;

	Random &random = Random::instance();

  double random_number = random.rand(); // between 0 and 1

	double E_sec = interpolate(random_number, dCS_interp, tabE_sec); // tabdCS have values between 0 and 1 (it's been normalised)

  double w = 1; // weigth of the particles



	if (haveElectrons && secondary == 0) {

    double f = (E_sec + mec2) / E;

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(11, (E_sec + mec2) / (1 + z) * eV / 1e-9, pos, w, interactionTag); // correct for the redshift
    	candidate->current.setEnergy((E - E_sec - mec2) / (1 + z) * eV / 1e-9 ); // E_sec is kinetic !! ; Em : rest energy
		}
	}

	else if (haveElectrons && secondary == 1) {

    double f = (E_sec + mec2) / E;

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			w = 1. / pow(1 - f, thinning); // w : weight of the secondary, default is 1
			candidate->addSecondary(-11, (E_sec + mec2) / (1 + z) * eV / 1e-9, pos, w, interactionTag);
      candidate->current.setEnergy((E - E_sec - mec2) / (1 + z) * eV / 1e-9 ); // E_sec is kinetic !! ; Em : rest energy
		}
	}

	else if (secondary == 2) {

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(1000010010, (E_sec + mpc2) / (1 + z) * eV / 1e-9, pos, w, interactionTag);
    candidate->current.setEnergy((E - E_sec - mpc2) / (1 + z) * eV / 1e-9); // E_sec is kinetic !! ; Em : rest energy
	}

	else if (secondary == 3) {

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(-1000010010, (E_sec + mpc2) / (1 + z) * eV / 1e-9, pos, w, interactionTag);
  	candidate->current.setEnergy((E - E_sec - mpc2) / (1 + z) * eV / 1e-9); // E_sec is kinetic !! ; Em : rest energy
	}

	else if (havePhotons && secondary == 4)  {

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(22, E_sec / (1 + z) * eV / 1e-9, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec) / (1 + z) * eV / 1e-9 ); // E_sec is kinetic !! ; Em : rest energy
	}


	// Update the primary particle energy.
	// This is done after adding the secondaries to correctly set the secondaries parent

  tabE_sec.clear();
	tabdCS.clear();

}

void HI::process(Candidate *candidate) const {

  //std::cout << "0!" << std::endl;


  int id = candidate->current.getId();

  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) / eV * 1e-9;

  if ((abs(id) != 1000010010) && (abs(id) != 1000020040) && (abs(id) != 1000060120) && (abs(id) != 1000130260) && (abs(id) != 1000260520 )) // we don't want particles that are not in the tables (H, He, C, Al, Fe)
		return;


  Vector3d position = candidate->current.getPosition();

  std::vector<double> array = setInteraction(candidate); // return 2 values, CS_tot in 0, and secondary in 1

  //std::cout << "10!" << std::endl;
  //std::cout << array[0] << " ";
  //std::cout << array[1] << " ";

  double inverse_mfp = density * array[0] * (barn * 1e-3)  ; //CS_tot computed before, put it in SI (from mb)

  std::cout << inverse_mfp << "  " << array[0] << "  " << density << std::endl;

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
		//std::cout << step << "  " << randDistance << std::endl;
		performInteraction(candidate, array[1]);
		step -= randDistance;
	} while (step > 0.);

  //std::cout << "OUT ! " << std::endl;
}




void HI::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string HI::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa



