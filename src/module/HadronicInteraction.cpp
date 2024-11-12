#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/massDistribution/ConstantDensity.h"

#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>


#include <iostream> //



namespace crpropa {

static const double mec2 = mass_electron * c_squared / GeV; // in GeV
static const double mpc2 = mass_proton * c_squared / GeV;
static const double mnc2 = mass_neutron * c_squared / GeV;



HadronicInteraction::HadronicInteraction( double H_density, double He_density, bool haveElectrons, bool havePhotons, bool haveNeutrinos, double thinning, double limit) {
  setGasDensity(H_density, He_density); // this one will probably disappear in order to be set from the simulation
	setHaveElectrons(haveElectrons);
	setHavePhotons(havePhotons);
	setHaveNeutrinos(haveNeutrinos);
	setLimit(limit);
	setThinning(thinning); // functions run to initialize the HadronicInteraction class
  setTables();
}


void HadronicInteraction::setGasDensity(double H_density, double He_density) {
	this-> H_density = H_density;
  this-> He_density = He_density;
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

void HadronicInteraction::setTables(){


  std::string name = "./filenames.txt";
  Filenames.clear();

	std::ifstream infile(name.c_str());

   part_numbering[1000010010] = "_p_";
   part_numbering[1000020040] = "_He_";
   part_numbering[1000060120] = "_C_";
   part_numbering[1000130260] = "_Al_";
   part_numbering[1000260520] = "_Fe_";

  if (!infile.good()){ // get the name of the data files
    throw std::runtime_error("HadronicInteraction: could not open file " + name);
    }

    while (infile.good()) {
      if (infile.peek() != '#') {
        std::string a;
        infile >> a;
        if (infile) {
          Filenames.push_back(a);
        }
      }
      infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n'); //to ignore the rest of the line
    }
  infile.close();


  for (std::string filename : Filenames) {   // collect the data in 3 dicts

    std::string path = "./Tables_HI/" + filename;

    if (filename.find("cum") != std::string::npos) {
        std::string reduced_name = filename;
        dict_cumDCS[reduced_name.erase(0, 4)] = std::vector<std::vector<double>>{};  // Erases k characters after idx, ie we get the same name as when it's not the cumulative cross sections but the real one
        dictE_sec[reduced_name] = std::vector<double>{};
        dictE_prim[reduced_name] = std::vector<double>{};

        std::ifstream infile(path.c_str());

        if (!infile.good()){
          throw std::runtime_error("HadronicInteraction: could not open file " + filename);
          }


        bool firstLine = true;
        double value;
        double firstValue;

        while (infile.good()) {
          if (infile.peek() != '#') {  // Ignore lines starting with '#'
              std::string line;
              if (std::getline(infile, line)) {
                  std::istringstream stream(line);


                  if (firstLine) {
                      stream >> firstValue;
                      while (stream >> value){
                        dictE_sec[reduced_name].push_back(value);
                      }

                      firstLine = false;
                  }

                  else{
                    stream >> firstValue;
                    dictE_prim[reduced_name].push_back(firstValue); // the list of Eprimary for cumulative data is smaller, because I separeted the low energy part from the rest (the Esec were different)
                    // Store the rest of the line into a row vector
                    std::vector<double> row;
                    while (stream >> value) {
                      row.push_back(value );  // store each remaining value ! Careful, it's mbarn
                    }

                    if (!row.empty()) {   // if row is not empty
                      dict_cumDCS[reduced_name].push_back(row); // each row correspond to a given primary energy, what varies are the sencondary energies that are being considered
                  }
                  }
              }
          }
          infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Skip to the next line
      }

    } else {
        dict[filename] = std::vector<std::vector<double>>{};
        dict[filename].resize(2); // create two empty vectors, so that the structure of the map is set

        //std::ifstream infile(("./Tables_HI/" + filename).c_str());
        std::ifstream infile(path.c_str());

        if (!infile.good()){
          throw std::runtime_error("HadronicInteraction: could not open file " + filename);
          }

          while (infile.good()) {
            if (infile.peek() != '#') {
              double a, b;
              infile >> a >> b;
              if (infile) {
                dict[filename][0].push_back(a); // a is the primary energy
                dict[filename][1].push_back(b); // b is the cross section
              }
            }
            infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n'); //to ignore the rest of the line
          }
        }
    infile.close();
  }

}



std::vector<double> HadronicInteraction::setInteraction(Candidate *candidate) const { // given a particle, get the primary total energies back, and use it (and the cross-sections) to decide what secondary will be produced (provided an interaction occurs). Then we need to decide whether the interaction occurs or not (given the total cross-section computed here, in process), and then perform the interaction (in performInteraction)

  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) / GeV ; // convert energy from J to GeV

  double CS_tot = 0;
  double secondary = 0;
  int index = 0;

  std::vector<double> tabCS_E; // table of cross sections (for different processes) at a given energy : used to decide which process we perform
  std::string chain = part_numbering.at(id); // can do any CR for which we have tables //  use .at() instead of just [] because we are in const functions, and part_numbering was defined outised this function
  std::vector<int> indexes;

  std::vector<std::string> tab_filenames;


  for (const auto& name : Filenames) { // get the files we are interested in (in that case the ones for which the primary is a proton)
    if (name.find(chain) != std::string::npos && name.find("04L") == std::string::npos &&  name.find("cum") == std::string::npos) { // we only want the interaction with protons for now, and the 04L are there only for cumulative tables, we don't want them now
        tab_filenames.push_back(name); // because it is an iterator
        indexes.push_back(index);
      }
    index += 1;
  }

  for (size_t i = 0; i < tab_filenames.size(); ++i) {

    tabCS_E.push_back(interpolate(E, dict.at(tab_filenames[i]).at(0), dict.at(tab_filenames[i]).at(1))); // get the CS at the energy of the primary interacting particle ; interpolate(x, tabX, tabY) ; we interpolate on the primary energy (x) and the cross section associated to a given process
    CS_tot += tabCS_E[i]; // -1 is not allowed in C++

  }

  Random &random = Random::instance();
  double random_number = random.rand() * CS_tot;
  double summed_CS = tabCS_E[0];

  while (summed_CS <= random_number){
    secondary += 1;
    summed_CS += tabCS_E[secondary];

  }
  return std::vector<double> {CS_tot, indexes[secondary]}; // indexes[secondary] is the index corresponding to the filename responsible for the interaction , ie we know the secondary being produced
}


void HadronicInteraction::performInteraction(Candidate *candidate, double secondary) const {

  std::string filename;

  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) / GeV;

  int filename_index = static_cast<int>(secondary); // get the int back so that we can acces the corresponding filename


  if (E > dict.at(Filenames[filename_index]).at(0).at(0) && E < 10){ // 10 GeV is the first energy of the primary for the non L files (ie the not low energy files) ; need the dict with the low E part
    std::string temp_filename = Filenames[filename_index]; // Make a non-const copy
    temp_filename.erase(temp_filename.size() - 4);         // Erase the last 4 characters
    filename = temp_filename + "L.txt";
  }
  else{
    filename = Filenames[filename_index]; // don't take the low energy file
  }


  std::vector<double> dCS_interp;


  const auto& data_row = dict_cumDCS.at(filename).at(0);

  for (size_t i = 0; i < data_row.size(); ++i) {
        std::vector<double> result;
        std::transform(dict_cumDCS.at(filename).begin(), dict_cumDCS.at(filename).end(), std::back_inserter(result),[i](const std::vector<double>& subvec) {return subvec[i];}); // extract ith element of each subvector of dict_cumDCS[filename] (the table containing the cumulative differential cross sections) // should be ok, just copied from previous code

        dCS_interp.push_back(interpolate(E, dictE_prim.at(filename), result)) ; // interpolate a cumulative differential cross section value, at a given secondary energy (we loop on those energy, via i), on primary energy, and differential cross section for different primary energies # since it's cumulative CS and that we add values at increasinigly high Esec, dCS_interp is ordered
    }


	Random &random = Random::instance();

  double random_number = random.rand(dCS_interp.back() - dCS_interp[0]) + dCS_interp[0]; // random.rand is written to generate random numbers between 0 and x ; so we generate a random number betwen 0 and x-min, and then add min to this random number ; back : last element of the vector


	double E_sec = interpolate(random_number, dCS_interp, dictE_sec.at(filename)); //  have values between 0 and 1 (it's been normalised) // at : because we try to access a constant variable, and either use it in loops or in other non constant functions, which could indeed modify this variable. Using dict.at(filename) instead means that we access this variable in read-only mode


  double w = 1; // weigth of the particles



	if (haveElectrons && filename.substr(0, 3) == "el_") { // ie the secondary is an electron

    double f = (E_sec + mec2) / E;

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(11, (E_sec + mec2) / (1 + z) * GeV , pos, w, interactionTag); // correct for the redshift
    	candidate->current.setEnergy((E - E_sec - mec2) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy  ; should also substract the energy even if haveElectrons == False
		}
	}

	else if (haveElectrons && filename.find("pos_") != std::string::npos) { // positron

    double f = (E_sec + mec2) / E;

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		if (random.rand() < pow(1 - f, thinning)) {
			w = 1. / pow(1 - f, thinning); // w : weight of the secondary, default is 1
			candidate->addSecondary(-11, (E_sec + mec2) / (1 + z) * GeV, pos, w, interactionTag);
      candidate->current.setEnergy((E - E_sec - mec2) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy
		}
	}

	else if (filename.find("prot_") != std::string::npos) { // proton

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(1000010010, (E_sec + mpc2) / (1 + z) * GeV, pos, w, interactionTag);
    candidate->current.setEnergy((E - E_sec - mpc2) / (1 + z) * GeV); // E_sec is kinetic !! ; Em : rest energy
	}

	else if (filename.find("aprot_") != std::string::npos) { // anti proton

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(-1000010010, (E_sec + mpc2) / (1 + z) * GeV, pos, w, interactionTag);
  	candidate->current.setEnergy((E - E_sec - mpc2) / (1 + z) * GeV); // E_sec is kinetic !! ; Em : rest energy
	}

	else if (havePhotons && filename.find("gam_") != std::string::npos)  { // gamma

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(22, E_sec / (1 + z) * GeV, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (haveNeutrinos && filename.find("nu_el") != std::string::npos)  { //  electronic neutrino, no mass energy for neutrinos ? We probably don't even care

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(12, E_sec / (1 + z) * GeV, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (haveNeutrinos && filename.find("anu_el") != std::string::npos)  { //  electronic antineutrino

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(-12, E_sec / (1 + z) * GeV, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (haveNeutrinos && filename.find("nu_mu") != std::string::npos)  { //  muonic neutrino

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(14, E_sec / (1 + z) * GeV, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (haveNeutrinos && filename.find("anu_mu") != std::string::npos)  { //  muonic antineutrino

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(-14, E_sec / (1 + z) * GeV, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (filename.substr(0, 2) == "n_")  { // neutron

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(2112, (E_sec + mnc2) / (1 + z) * GeV, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec - mnc2) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (filename.find("an_") != std::string::npos)  { //  anti neutron

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(-2112, (E_sec + mnc2) / (1 + z) * GeV, pos, w, interactionTag); // photons # Em = 0
  	candidate->current.setEnergy((E - E_sec - mnc2) / (1 + z) * GeV ); // E_sec is kinetic !! ; Em : rest energy  // got all our particles for now, need to add the mass energy of each particle
	}

	dCS_interp.clear();

}

void HadronicInteraction::process(Candidate *candidate) const {


  int id = candidate->current.getId();

  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) / GeV ;

  if ((abs(id) != 1000010010) && (abs(id) != 1000020040) && (abs(id) != 1000060120) && (abs(id) != 1000130260) && (abs(id) != 1000260520 )) // we don't want particles that are not in the tables (H, He, C, Al, Fe)
		return;


  Vector3d position = candidate->current.getPosition();

  std::vector<double> array = setInteraction(candidate); // return 2 values, CS_tot in 0, and secondary in 1

  double inverse_mfp = (H_density + He_density) * array[0] * (barn * 1e-3)  ; //CS_tot computed before, put it in SI (from mb) // this one is not updated if we do several interactions in one step (the energy of the primary changes)


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
		performInteraction(candidate, array[1]);
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
