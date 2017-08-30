#include<vector>
#include<iostream>
#include"TempSwimSpeedOptim.hpp"

int main(void) {
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::exp_ratio_fitness;
	using predator_fitness = tempss::difference_fitness;
	using this_system = tempss::stepdrdm_optimizer_system<prey_fitness, predator_fitness>;
	using this_state = typename this_system::state;
	using freq_state = tempss::freq_state;
	std::vector<double> V{ 1.000, 1.298, 1.467, 1.651, 1.848, 2.061, 2.290, 2.536, 2.798, 3.077, 3.375};
	std::vector<double> U{ 1.000, 1.091, 1.136, 1.181, 1.227, 1.283, 1.318, 1.364, 1.409, 1.454, 1.500};
	std::vector<double> C{ 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100};
	double a = 0.5;
	double b = 2.0;
	double h = 2.0;
	double k = 0.1;
	double c = 0.0001;
	double bc = 0.1;
	double d = 0.01;
	double e = 0.3;

	this_system System(predation(a, b, h), prey_fitness(1.0, 0.0), predator_fitness(), V.begin(), V.end(), U.begin(), U.end(), C.begin(), C.end(), d, e, k);

	this_state vLower;
	this_state vUpper;
	std::tie(vLower, vUpper) = System();

	freq_state Lower = System.get_freq_state(vLower);
	freq_state Upper = System.get_freq_state(vUpper);

	system("pause");
	return 0;
}
