#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>



namespace crpropa {

static const double mec2 = mass_electron * c_squared ; // in J
static const double mpc2 = mass_proton * c_squared ;
static const double mnc2 = mass_neutron * c_squared ;



HadronicInteraction::HadronicInteraction( double H_density, double He_density, bool haveElectrons, bool havePhotons, bool haveNeutrinos, bool haveantiNucleons, double thinning, double limit) {
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

void HadronicInteraction::setHaveantiNucleons(bool haveantiNucleons) {
	this->haveantiNucleons = haveantiNucleons;
}

void HadronicInteraction::setLimit(double limit) {
	this->limit = limit;
}

void HadronicInteraction::setThinning(double thinning) {
	this->thinning = thinning;
}

void HadronicInteraction::setTables(){

  data_map[1000010010] = std::map<int, std::vector<std::vector<std::vector<double>>>>(); // first layer of the map, based on the CR identity
  data_map[1000020040] = std::map<int, std::vector<std::vector<std::vector<double>>>>();
  data_map[1000060120] = std::map<int, std::vector<std::vector<std::vector<double>>>>();
  data_map[1000130260] = std::map<int, std::vector<std::vector<std::vector<double>>>>();
  data_map[1000260520] = std::map<int, std::vector<std::vector<std::vector<double>>>>();
  data_map[1000010012] = std::map<int, std::vector<std::vector<std::vector<double>>>>(); // + 2 at the end : to say that the target is Helium
  data_map[1000020042] = std::map<int, std::vector<std::vector<std::vector<double>>>>();


  std::string name = "./filenames.txt";
  Filenames.clear();

	std::ifstream infile(name.c_str());


  // get the name of the data files
  if (!infile.good()){
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


  int id; // identity of the CR
  int secondary;

  for (std::string filename : Filenames) {   // collect the data in a map of map of tuple of vectors (and then I can only use ints) and keep track of the primary and secondary

    std::string path = "./Tables_HI_extra/" + filename;


    // initialize the maps
    for (const auto& [key_CR, value_CR] : CR_numbering) {
        if (filename.find(key_CR) != std::string::npos) {
          if (filename.find("p04")!= std::string::npos) {
            id = value_CR;
          } else {
              id = value_CR + 2; // ie the target is helium
            }
          break; // avoid continuing in the loop for nothing
        }
      }
      std::string starter;
      if (filename.find("cum") != std::string::npos ) {
        starter = filename.substr(4, 6); // 6 : length of the string
      }
      else {
        starter = filename.substr(0, 6);
      }
      for (const auto& [key, value] : secondaries_numbering) {
        if (starter.find(key) == 0) { // otherwise I had some issues with for example "n_" and "an_" : the second contains the first, hence find is not good enough, we also need the location
          if (data_map.find(id) == data_map.end() || data_map[id].find(value) == data_map[id].end()) { // check if we already initialized it

            data_map[id][value] = std::vector<std::vector<std::vector<double>>>();
            data_map[id][value].resize(5);

            data_map[id][value][0].resize(2);
            data_map[id][value][1].resize(2);
            data_map[id][value][3].resize(2);

          }
          secondary = value;
          break;
        }
      }
    int to_fill1 = 1;
    int to_fill2 = 2;

    if (filename.find("04L") != std::string::npos){
      to_fill1 = 3;   // ie in that case we fill in the two last vectors, where we put the (Eprimary Esecondary) we have for this file, and the differential cross section
      to_fill2 = 4;
    }

    //std::cout << id << "  " << secondary << "  " << to_fill2 << "  " << filename << std::endl;

    if (filename.find("cum") != std::string::npos ) {


      std::ifstream infile(path.c_str());

      if (!infile.good()){
        throw std::runtime_error("HadronicInteraction: could not open file " + filename);
        }

      bool firstLine = true;
      double data;
      double firstValue;

        while (infile.good()) {
          if (infile.peek() != '#') {  // Ignore lines starting with '#'
              std::string line;
              if (std::getline(infile, line)) {
                  std::istringstream stream(line);

                  if (firstLine) {
                      stream >> firstValue;
                      while (stream >> data){
                        data_map[id][secondary][to_fill1][1].push_back(data * GeV);  // put it in Joules // st::get : needed because this is a tuple // energies of secondaries
                      }

                      firstLine = false;
                  }

                  else{
                    stream >> firstValue;
                    data_map[id][secondary][to_fill1][0].push_back(firstValue * GeV); // the list of Eprimary for cumulative data is smaller, because I separeted the low energy part from the rest (the Esec were different)
                    // Store the rest of the line into a row vector
                    std::vector<double> row;
                    while (stream >> data) {
                      row.push_back(data);  // store each remaining value ! No unit as I normalised it (ie the tables are normalised between 0 and 1)
                    }

                    if (!row.empty()) {   // if row is not empty
                      data_map[id][secondary][to_fill2].push_back(row); // each row correspond to a given primary energy, what varies are the secondary energies that are being considered
                      const auto& lol = data_map[id][secondary][to_fill2][0];
                    }
                  }
              }
          }
          infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Skip to the next line
      }

    } else {

        std::ifstream infile(path.c_str());

        if (!infile.good()){
          throw std::runtime_error("HadronicInteraction: could not open file " + filename);
          }

          while (infile.good()) {
            if (infile.peek() != '#') {
              double a, b;
              infile >> a >> b;
              if (infile) {
                data_map[id][secondary][0][0].push_back(a * GeV); // a is the primary energy, put it in J
                data_map[id][secondary][0][1].push_back(b * (barn * 1e-3)); // b is the cross section ; again, from mbarns to SI (m^2)
              }
            }
            infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n'); //to ignore the rest of the line
          }
        }

    infile.close();
  }
std::cout << "Data loaded" << std::endl;
}



std::vector<double> HadronicInteraction::setInteraction(Candidate *candidate) const { // given a particle, gets the primary total energies back, and uses it (and the cross-sections) to decide what secondary will be produced (provided an interaction occurs). Then we need to decide whether the interaction occurs or not (given the total cross-section computed here, in process), and then perform the interaction (in performInteraction)


  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) ;

  double CS_H = 0;
  double CS_He = 0;
  int secondary = 0;

  std::vector<double> tabCS_E; // table of cross sections (for different processes) at a given energy : used to decide which process we perform

  for (size_t i = 0; i < secondaries.size(); ++i) {

    tabCS_E.push_back(interpolate(E, data_map.at(id).at(secondaries.at(i)).at(0).at(0), data_map.at(id).at(secondaries.at(i)).at(0).at(1))) ; // get the CS at the energy of the primary interacting particle ; interpolate(x, tabX, tabY) ; we interpolate on the primary energy (x) and the cross section associated to a given process

    CS_H += tabCS_E.back(); // -1 is not allowed in C++

    tabCS_E.back() *= H_density;

  }

  if ((id == 1000010010 ) || (id == 1000020040)) {

    for (size_t i = 0; i < secondaries.size(); ++i) { // do the same but with id+2, to account for He as a target

      tabCS_E.push_back(interpolate(E, data_map.at(id + 2).at(secondaries.at(i)).at(0).at(0), data_map.at(id + 2).at(secondaries.at(i)).at(0).at(1))); // get the CS at the energy of the primary interacting particle ; interpolate(x, tabX, tabY) ; we interpolate on the primary energy (x) and the cross section associated to a given process

      CS_He += tabCS_E.back(); // -1 is not allowed in C++

      tabCS_E.back() *= He_density;
      }
    }

  
  Random &random = Random::instance();
  double random_number = random.rand() * (CS_H * H_density + CS_He * He_density); // need to account for densities (if there is no He, we can’t interact with it)
  double summed_CS = tabCS_E[0];

  while (summed_CS <= random_number){
    secondary += 1;
    summed_CS += tabCS_E[secondary];

  }
  return std::vector<double> {CS_H, CS_He, secondary}; //  ie we keep track of the secondary being produced

}


void HadronicInteraction::performInteraction(Candidate *candidate, double sec) const {

  int id = candidate->current.getId();
  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) ;
  int to_read = 0;

  int secondary = static_cast<int>(sec); // get the int back so that we can acces the corresponding filename

  if (((id == 1000010010) || (id == 1000020040))  && (E < 1.602176487e-9)){ // 10 GeV = 1.602176487e-9 is the first energy of the primary for the non L files (ie the not low energy files) ; the low energy extension is only for proton and He
    to_read += 2;
    }

  if (secondary >=  secondaries.size()){
    secondary -= secondaries.size(); // ie we have a CR interacting with He
    id += 2; // to read the correct tables
  }

  secondary = secondaries.at(secondary); // get the code of the secondaries being produced

  std::vector<double> dCS_interp;
  std::vector<double> result;

  const auto& data_row = data_map.at(id).at(secondary).at(2 + to_read).at(0); // check the indexes maybe


  for (size_t i = 0; i < data_row.size(); ++i) {

        //std::transform(dict_cumDCS.at(filename).begin(), dict_cumDCS.at(filename).end(), std::back_inserter(result),[i](const std::vector<double>& subvec) {return subvec[i];}); // extract ith element of each subvector of dict_cumDCS[filename] (the table containing the cumulative differential cross sections) // should be ok, just copied from previous code
        //std::cout << "about to enter the second for" << std::endl;
        for (auto it = data_map.at(id).at(secondary).at(2 + to_read).begin(); it != data_map.at(id).at(secondary).at(2 + to_read).end(); ++it) {
          //std::cout << "in the second for" << std::endl;
          result.push_back((*it)[i]);
        }

        dCS_interp.push_back(interpolate(E, data_map.at(id).at(secondary).at(1 + to_read).at(0), result)) ; // interpolate a cumulative differential cross section value, at a given secondary energy (we loop on those energy, via i), on primary energy, and differential cross section for different primary energies # since it's cumulative CS and that we add values at increasinigly high Esec, dCS_interp is ordered
        result.clear();
    }


	Random &random = Random::instance();

  double random_number = random.rand(dCS_interp.back() - dCS_interp[0]) + dCS_interp[0]; // random.rand is written to generate random numbers between 0 and x ; so we generate a random number betwen 0 and x-min, and then add min to this random number ; back : last element of the vector // maybe this could be simplified

	double E_sec = interpolate(random_number, dCS_interp, data_map.at(id).at(secondary).at(1 + to_read).at(1)); //  have values between 0 and 1 (it's been normalised) // at : because we try to access a constant variable, and either use it in loops or in other non constant functions, which could indeed modify this variable. Using dict.at(filename) instead means that we access this variable in read-only mode



  double w = 1; // weigth of the particles



	if (secondary == 11) { // ie the secondary is an electron

    candidate->current.setEnergy((E - E_sec - mec2) / (1 + z) ); // E_sec is kinetic !! ; Em : rest energy  ; should also substract the energy even if haveElectrons == False

    if (haveElectrons){
      double f = (E_sec + mec2) / E;

      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      if (random.rand() < pow(1 - f, thinning)) {
        w = 1. / pow(1 - f, thinning);
        candidate->addSecondary(11, (E_sec + mec2) / (1 + z) , pos, w, interactionTag); // correct for the redshift
      }
		}
	}

	else if (secondary == -11) { // positron

    candidate->current.setEnergy((E - E_sec - mec2) / (1 + z) ); // E_sec is kinetic !! ; Em : rest energy

    if (haveElectrons){
      double f = (E_sec + mec2) / E;

      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      if (random.rand() < pow(1 - f, thinning)) {
        w = 1. / pow(1 - f, thinning); // w : weight of the secondary, default is 1
        candidate->addSecondary(-11, (E_sec + mec2) / (1 + z) , pos, w, interactionTag);
      }
		}
	}

	else if (secondary == 1000010010) { // proton

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(secondary, (E_sec + mpc2) / (1 + z), pos, w, interactionTag);
    candidate->current.setEnergy((E - E_sec - mpc2) / (1 + z)); // E_sec is kinetic !! ; Em : rest energy
	}

	else if (secondary == -1000010010) { // anti proton

    candidate->current.setEnergy((E - E_sec - mpc2) / (1 + z) ); // E_sec is kinetic !! ; Em : rest energy

    if (haveantiNucleons){
      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      candidate->addSecondary(secondary, (E_sec + mpc2) / (1 + z), pos, w, interactionTag);
    }
	}

	else if (secondary == 22)  { // gamma
    if (havePhotons){
      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      candidate->addSecondary(secondary, E_sec / (1 + z) , pos, w, interactionTag); // photons # Em = 0
    }
  	candidate->current.setEnergy((E - E_sec) / (1 + z)  ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (secondary == 12)  { //  electronic neutrino, no mass energy for neutrinos ? We probably don't even care
    candidate->current.setEnergy((E - E_sec) / (1 + z)  ); // E_sec is kinetic !! ; Em : rest energy

    if (haveNeutrinos){
      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      candidate->addSecondary(secondary, E_sec / (1 + z), pos, w, interactionTag);
    }
	}

  else if (secondary == -12)  { //  electronic neutrino, no mass energy for neutrinos ? We probably don't even care
    candidate->current.setEnergy((E - E_sec) / (1 + z)  ); // E_sec is kinetic !! ; Em : rest energy

    if (haveNeutrinos){
      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      candidate->addSecondary(secondary, E_sec / (1 + z), pos, w, interactionTag);
    }
	}

  else if (secondary == 14)  { //  electronic neutrino, no mass energy for neutrinos ? We probably don't even care
    candidate->current.setEnergy((E - E_sec) / (1 + z)  ); // E_sec is kinetic !! ; Em : rest energy

    if (haveNeutrinos){
      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      candidate->addSecondary(secondary, E_sec / (1 + z), pos, w, interactionTag);
    }
	}

  else if (secondary == -14)  { //  electronic neutrino, no mass energy for neutrinos ? We probably don't even care
    candidate->current.setEnergy((E - E_sec) / (1 + z)  ); // E_sec is kinetic !! ; Em : rest energy

    if (haveNeutrinos){
      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      candidate->addSecondary(secondary, E_sec / (1 + z), pos, w, interactionTag);
    }
	}

  else if (secondary == 100000010)  { // neutron

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->addSecondary(secondary, (E_sec + mnc2) / (1 + z) , pos, w, interactionTag);
  	candidate->current.setEnergy((E - E_sec - mnc2) / (1 + z)  ); // E_sec is kinetic !! ; Em : rest energy
	}

  else if (secondary == -100000010)  { //  anti neutron

    candidate->current.setEnergy((E - E_sec - mnc2) / (1 + z)  ); // E_sec is kinetic !! ; Em : rest energy  // got all our particles for now

    if (haveantiNucleons){
      Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
      candidate->addSecondary(secondary, (E_sec + mnc2) / (1 + z) , pos, w, interactionTag);
    }

	}
}

void HadronicInteraction::process(Candidate *candidate) const {


  int id = candidate->current.getId();

  double z = candidate->getRedshift();
  double E = candidate->current.getEnergy() * (1+z) ;


  if ((id != 1000010010) && (id != 1000020040) && (id != 1000060120) && (id != 1000130260) && (id != 1000260520 )) // we don't want particles that are not in the tables (H, He, C, Al, Fe)
		return;


  if ((E > data_map.at(id).at(11).at(0).at(0).back()) || (E < data_map.at(id).at(11).at(0).at(0).at(0) )) { // all the limits are secondaries independant, apart from exotic ones (ad for example)
    return;
  }

  Vector3d position = candidate->current.getPosition();

  std::vector<double> array = setInteraction(candidate); // return 3 values, CS with p and He in 0 and 1, and the code fo the secondary in 2

  double inverse_mfp = H_density * array[0] + He_density * array[1]; //CS_tot computed before, // this one is not updated if we do several interactions in one step (the energy of the primary changes)


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
		performInteraction(candidate, array[2]);
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
