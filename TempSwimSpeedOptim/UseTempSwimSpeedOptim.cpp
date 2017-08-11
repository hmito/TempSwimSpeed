#include<vector>
#include<iostream>
#include"TempSwimSpeedOptim.hpp"

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

	auto x1 = System.optimize_by_hill_climbing(1000,100000);
	double W1 = System.get_prey_fitness(x1);

	boost::format F("%.3f\t%.3f\t%.3f\t%.3f\n");

	std::cout << "===x1===" << std::endl;
	for (unsigned int i = 0; i < U.size(); ++i) {
		double f = x1[i];
		double r = System.begin()[i].prey_reward(f);
		double m = System.begin()[i].prey_mortality(f);
		double drdm = r / (m+std::numeric_limits<double>::min());

		std::cout << F%f%r%m%drdm << std::endl;
	}
	std::cout << std::endl;

	auto x2 = System.optimize_by_small_step();
	double W2 = System.get_prey_fitness(x2);

	std::cout << "===x2===" << std::endl;
	for (unsigned int i = 0; i < U.size(); ++i) {
		double f = x2[i];
		double r = System.begin()[i].prey_reward(f);
		double m = System.begin()[i].prey_mortality(f);
		double drdm = r / (m + std::numeric_limits<double>::min());

		std::cout << F%f%r%m%drdm << std::endl;
	}
	std::cout << std::endl;

	auto x3 = System.optimize_by_bisect();
	double W3l = System.get_prey_fitness(x3.first);
	double W3u = System.get_prey_fitness(x3.second);

	std::cout << "===x3l===" << std::endl;
	for (unsigned int i = 0; i < U.size(); ++i) {
		double f = x3.first[i];
		double r = System.begin()[i].prey_reward(f);
		double m = System.begin()[i].prey_mortality(f);
		double drdm = r / (m + std::numeric_limits<double>::min());

		std::cout << F%f%r%m%drdm << std::endl;
	}

	std::cout << "===x3u===" << std::endl;
	for (unsigned int i = 0; i < U.size(); ++i) {
		double f = x3.second[i];
		double r = System.begin()[i].prey_reward(f);
		double m = System.begin()[i].prey_mortality(f);
		double drdm = r / (m + std::numeric_limits<double>::min());

		std::cout << F%f%r%m%drdm << std::endl;
	}
	std::cout << std::endl;

	system("pause");
	return 0;
}
