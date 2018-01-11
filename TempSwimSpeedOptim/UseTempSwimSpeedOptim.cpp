#include<vector>
#include<iostream>
#include"TempSwimSpeedOptim.hpp"

int main(void) {
	std::vector<double> V{ 1.93473690388997, 1.80865828381746, 1.69561928549564, 1.60332332985438, 1.53806023374436, 1.50427756931309, 1.50427756931309, 1.53806023374436, 1.60332332985438, 1.69561928549564, 1.80865828381746, 1.93473690388997, 2.06526309611003, 2.19134171618254, 2.30438071450436, 2.39667667014562, 2.46193976625564, 2.49572243068691, 2.49572243068691, 2.46193976625564, 2.39667667014562, 2.30438071450436, 2.19134171618254, 2.06526309611003 };
	std::vector<double> U{ 0.654984994781574, 0.557090350616535, 0.506416353969642, 0.506416353969642, 0.557090350616535, 0.654984994781574, 0.793428928243459, 0.962987425726183, 1.15210535583496, 1.34789464416504, 1.53701257427382, 1.70657107175654, 1.84501500521843, 1.94290964938346, 1.99358364603036, 1.99358364603036, 1.94290964938346, 1.84501500521843, 1.70657107175654, 1.53701257427382, 1.34789464416504, 1.15210535583496, 0.962987425726183, 0.793428928243459 };
	std::vector<double> K{ 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02 };
	std::vector<double> C{ 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015 };
	std::vector<double> L{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
//	double a = 0.5;
	double d = 0.001;
	double e = 0.25;
	double b = 0.0;
	double h = 2.0;
//	double k = 0.1;
//	double c = 0.0001;
//	double base_c = 0.1;
	double cb = 1e-4;
	double cf = 1e-4;

	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::exp_ratio_fitness;
	using predator_fitness = tempss::difference_fitness;
	using this_system = tempss::stepdrdm_optimizer_system<prey_fitness, predator_fitness>;
	using this_state = typename this_system::state;
	using freq_state = tempss::freq_state;

	this_system System(predation(1.0, b, h), prey_fitness(1.0, 0.0), predator_fitness(), V.begin(), V.end(), U.begin(), U.end(), K.begin(), K.end(), C.begin(), C.end(), L.begin(), L.end(), d, e, cb, cf);

	this_state svLower;
	this_state svUpper;
	std::tie(svLower, svUpper) = System();

	freq_state vLower = System.get_freq_state(svLower);
	freq_state vUpper = System.get_freq_state(svUpper);

	double PreyWL = System.get_prey_fitness(svLower);
	double PreyWH = System.get_prey_fitness(svUpper);
	double PredatorWL = System.get_predator_fitness(svLower);
	double PredatorWH = System.get_predator_fitness(svUpper);

	system("pause");
	return 0;
}
