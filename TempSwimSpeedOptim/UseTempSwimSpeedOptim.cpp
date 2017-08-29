#include<vector>
#include<iostream>
#include"TempSwimSpeedOptim.hpp"

int main(void) {
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::ratio_fitness;
	using predator_cost = tempss::constant_predator_cost;
	using prey_reward = tempss::linear_prey_reward;
	using this_system = tempss::optimizer_system<predation, prey_fitness>;

	std::vector<double> V{ 1.00,1.298,1.467,1.651,1.848,2.061,2.290,2.536,2.798,3.077,3.375};
	std::vector<double> U{ 1.00,1.091,1.136,1.181,1.227,1.283,1.318,1.364,1.409,1.454,1.500};
	double a = 0.5;
	double b = 2.0;
	double h = 2.0;
	double k = 0.1;
	double mu_a = 0.0005;
	double mu_r = 0.0005;
	double c = 0.0001;
	double d = 0.01;

	this_system System(predation(a, b, h), prey_fitness(1.0, V.size()*mu_r), prey_reward(k), predator_cost(c), V.begin(), V.end(), U.begin(), U.end(), mu_a - mu_r, d);

	tempss::state vLower;
	tempss::state vUpper;
	std::tie(vLower, vUpper) = System.optimize_by_stepbisect();

	system("pause");
	return 0;
}
