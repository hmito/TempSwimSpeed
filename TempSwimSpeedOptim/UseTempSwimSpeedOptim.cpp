#include<vector>
#include<iostream>
#include"TempSwimSpeedOptim.hpp"

int main(void) {
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::exp_diff_fitness;
	using predator_cost = tempss::constant_predator_cost;
	using prey_reward = tempss::linear_prey_reward;
	using this_system = tempss::optimizer_system<predation, prey_fitness>;


	std::vector<double> U{ 0.10,0.40,0.55,0.70,0.80,0.85, 0.90, 0.95 };
	std::vector<double> V{ 0.05,0.40,0.60,0.90,1.20,1.50,1.70, 1.75 };
	double alpha = 1.2;
	double n = 2.0;
	double h = 2.0;
	double beta = 4.0;
	double k = 2.0;
	double c = 1.0;
	double mu = 0.2;
	double d = 0.5;



	this_system System(predation(alpha, n, h), prey_fitness(beta), prey_reward(k), predator_cost(c), V.begin(), V.end(), U.begin(), U.end(), mu, d);

	auto Ans = System.optimize_by_stepbisect();

	tempss::state vLower;
	tempss::state vUpper;
	double Freq;
	double W;

	system("pause");
	return 0;
}
