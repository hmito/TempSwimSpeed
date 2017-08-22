//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include "TempSwimSpeedOptim/TempSwimSpeedOptim.hpp"

//probability foraging optimization
//V,U		Vector of speed of predators and prey at each time step
//following three parameters determine the predation rate: a*(v-u)^b / {1 + h*a*(v-u)^b} 
//	a		inverse of searching time	
//	b		non-linear influence of speed difference
//	h		average handling time for predation a prey
//following three parameters determine the prey traits
//	k		coefficient of foraging reward for prey (k*u is the reward)
//	mu_a	basic mortality rate of active prey (not include the mortality by predation)
//	mu_r	basic mortality rate of resting prey (mu_a >= mu_r)
//following two parameters determine the predator traits
//	d		relative density of predator/prey
//	c		predation cost for predators
// [[Rcpp::export]]
Rcpp::List tss_probforage_optimize(
	Rcpp::NumericVector V,
	Rcpp::NumericVector U,
	double a,
	double b,
	double h,
	double k,
	double mu_a,
	double mu_r,
	double c,
	double d
){
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::ratio_fitness;
	using predator_cost = tempss::constant_predator_cost;
	using prey_reward = tempss::linear_prey_reward;
	using this_system = tempss::optimizer_system<predation, prey_fitness>;

	this_system System(predation(a, b, h), prey_fitness(1.0, V.size()*mu_r), prey_reward(k), predator_cost(c), V.begin(), V.end(), U.begin(), U.end(), mu_a-mu_r, d);

	tempss::state vLower;
	tempss::state vUpper;
	std::tie(vLower,vUpper) = System();

	double PreyWL = System.get_prey_fitness(vLower);
	double PreyWH = System.get_prey_fitness(vUpper);
	double PredatorWL = System.get_predator_fitness(vLower);
	double PredatorWH = System.get_predator_fitness(vUpper);

	Rcpp::NumericVector PreyL(vLower.begin(), vLower.end());
	Rcpp::NumericVector PreyH(vUpper.begin(), vUpper.end());

	Rcpp::NumericVector PredatorL(vLower.size());
	Rcpp::NumericVector PredatorH(vLower.size());

	Rcpp::NumericVector rf(vLower.size());
	Rcpp::NumericVector pm1(vLower.size());
	Rcpp::NumericVector Thr(vLower.size());
	Rcpp::NumericVector drdm_f(vLower.size());
	Rcpp::NumericVector drdm_1(vLower.size());
	for(unsigned int i = 0; i < System.size(); ++i){
		const auto& TimeInfo = System.at(i);
		PredatorL[i] = TimeInfo.predator_strategy(PreyL[i]);
		PredatorH[i] = TimeInfo.predator_strategy(PreyH[i]);
		Thr[i] = TimeInfo.f_threshold();
		rf[i] = TimeInfo.prey_reward(0);
		pm1[i] = TimeInfo.prey_predation_mortality(1);
		drdm_f[i] = System.at(i).prey_drdm0();
		drdm_1[i] = System.at(i).prey_drdm1();
	}

	return Rcpp::List::create(
		Rcpp::Named("PreyL") = PreyL,
		Rcpp::Named("PreyH") = PreyH,
		Rcpp::Named("PreyWL") = PreyWL,
		Rcpp::Named("PreyWH") = PreyWH,
		Rcpp::Named("PredatorL") = PredatorL,
		Rcpp::Named("PredatorH") = PredatorH,
		Rcpp::Named("PredatorWL") = PredatorWL,
		Rcpp::Named("PredatorWH") = PredatorWH,
		Rcpp::Named("ThresholdPreyFreq") = Thr,
		Rcpp::Named("PreyReward") = rf,
		Rcpp::Named("PredationMortality_1") = pm1,
		Rcpp::Named("drdm_f") = drdm_f,
		Rcpp::Named("drdm_1") = drdm_1
	);
}
