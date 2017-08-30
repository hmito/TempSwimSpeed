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
//	mu_i	basic mortality rate of resting prey (mu_a >= mu_r)
//following two parameters determine the predator traits
//	M		Vector of predation cost for predators
//	d		relative density of predator/prey
// [[Rcpp::export]]
Rcpp::List tss_probforage_optimize(
	Rcpp::NumericVector V,
	Rcpp::NumericVector U,
	Rcpp::NumericVector C,
	double a,
	double b,
	double h,
	double k,
	double d,
	double e
){
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::exp_ratio_fitness;
	using predator_fitness = tempss::difference_fitness;
	using this_system = tempss::stepdrdm_optimizer_system<prey_fitness, predator_fitness>;
	using this_state = typename this_system::state;
	using freq_state = tempss::freq_state;

	this_system System(predation(a, b, h), prey_fitness(1.0, 0.0), predator_fitness(), V.begin(), V.end(), U.begin(), U.end(), C.begin(), C.end(), d, e, k);

	this_state svLower;
	this_state svUpper;
	std::tie(svLower, svUpper) = System();

	freq_state vLower = System.get_freq_state(svLower);
	freq_state vUpper = System.get_freq_state(svUpper);

	double PreyWL = System.get_prey_fitness(svLower);
	double PreyWH = System.get_prey_fitness(svUpper);
	double PredatorWL = System.get_predator_fitness(svLower);
	double PredatorWH = System.get_predator_fitness(svUpper);

	Rcpp::NumericVector PreyL(vLower.begin(), vLower.end());
	Rcpp::NumericVector PreyH(vUpper.begin(), vUpper.end());

	Rcpp::NumericVector PredatorL(vLower.size());
	Rcpp::NumericVector PredatorH(vLower.size());

	Rcpp::NumericVector rf(vLower.size());
	Rcpp::NumericVector Thr(vLower.size());
	Rcpp::NumericVector m0(vLower.size());
	Rcpp::NumericVector mF(vLower.size());
	Rcpp::NumericVector m1(vLower.size());
	for(unsigned int i = 0; i < System.size(); ++i){
		const auto& TimeInfo = System.at(i);
		PredatorL[i] = TimeInfo.predator_strategy(svLower[i]);
		PredatorH[i] = TimeInfo.predator_strategy(svUpper[i]);
		Thr[i] = TimeInfo.f_threshold();
		rf[i] = TimeInfo.prey_reward(0);
		m0[i] = TimeInfo.prey_mortality(0);
		mF[i] = TimeInfo.prey_mortality(1);
		m1[i] = TimeInfo.prey_mortality(2);
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
		Rcpp::Named("PreyMortality0") = m0,
		Rcpp::Named("PreyMortalityF") = mF,
		Rcpp::Named("PreyMortality1") = m1
	);
}
