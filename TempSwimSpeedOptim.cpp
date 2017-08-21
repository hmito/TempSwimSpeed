//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include "TempSwimSpeedOptim/TempSwimSpeedOptim.hpp"

//probability foraging optimization
//V,U	Speed of predators,prey
//alpha	inverse of searching time	
//n		non-linear influence of speed difference
//h		average handling time for predation a prey
//beta	fitness saturation for prey
//k		coefficient of foraging reward for prey
//c		predation cost for predators
//mu	basic mortality for prey
//d		relative density of predator/prey
// [[Rcpp::export]]
Rcpp::List tss_probforage_optimize(
	Rcpp::NumericVector V,
	Rcpp::NumericVector U,
	double alpha,
	double beta,
	double h,
	double lambda,
	double k,
	double c,
	double mu,
	double d
){
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::exp_diff_fitness;
	using predator_cost = tempss::constant_predator_cost;
	using prey_reward = tempss::linear_prey_reward;
	using this_system = tempss::optimizer_system<predation, prey_fitness>;

	this_system System(predation(alpha, beta, h), prey_fitness(lambda), prey_reward(k), predator_cost(c), V.begin(), V.end(), U.begin(), U.end(), mu, d);

	tempss::state vLower;
	tempss::state vUpper;
	double Freq;
	double W;
	std::tie(vLower,vUpper,Freq,W) = System.optimize_by_stepbisect();

	Rcpp::NumericVector PreyL(vLower.begin(), vLower.end());
	Rcpp::NumericVector PreyH(vUpper.begin(), vUpper.end());

	Rcpp::NumericVector PredatorL(vLower.size());
	Rcpp::NumericVector PredatorH(vLower.size());
	Rcpp::NumericVector Thr(vLower.size());
	Rcpp::NumericVector drdmL(vLower.size());
	Rcpp::NumericVector drdmH(vLower.size());
	for(unsigned int i = 0; i < System.size(); ++i){
		PredatorL[i] = System.at(i).predator_strategy(PreyL[i]);
		PredatorH[i] = System.at(i).predator_strategy(PreyH[i]);
		Thr[i] = System.at(i).f_threshold();
		drdmL[i] = System.at(i).prey_drdm0();
		drdmH[i] = System.at(i).prey_drdm1();
	}

	return Rcpp::List::create(
		Rcpp::Named("PreyL") = PreyL,
		Rcpp::Named("PreyH") = PreyH,
		Rcpp::Named("PredatorL") = PredatorL,
		Rcpp::Named("PredatorH") = PredatorH,
		Rcpp::Named("Thr") = Thr,
		Rcpp::Named("drdmL") = drdmL,
		Rcpp::Named("drdmH") = drdmH,
		Rcpp::Named("Freq") = Freq,
		Rcpp::Named("W") = W
	);
}
