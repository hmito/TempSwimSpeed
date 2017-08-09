#include"TempSwimSpeedOptim.hpp"
#include<vector>
#include<iostream>

int main(void) {
	using predation = tempss::search_hunt_forage_predation;
	using prey_fitness = tempss::exp_diff_fitness;
	using predator_cost = tempss::constant_predator_cost;
	using prey_reward = tempss::linear_prey_reward;
	using this_system = tempss::optimizer_system<predation, prey_fitness>;


	double s = 1;
	double d = 2;
	double h = 2;
	double alpha = 0.5;
	double k = 2.0;
	double c = 0.1;
	std::vector<double> U{ 0.10,0.40,0.55,0.70,0.80,0.85, 0.90, 0.95 };
	std::vector<double> V{ 0.05,0.40,0.60,0.90,1.20,1.50,1.70, 1.75 };
	double mu = 0.1;
	double rd = 0.1;

	this_system System(predation(s, d, h), prey_fitness(alpha), prey_reward(k), predator_cost(c), V.begin(), V.end(), U.begin(), U.end(), mu, rd);

//	auto x = System();
	auto x1 = System.optimize_by_hill_climbing();
	auto x2 = System.optimize_by_small_step();


	double W1 = System.get_prey_fitness(x1);

	double W2 = System.get_prey_fitness(x2);

	system("pause");
	return 0;
}