//TempSwimSpeed_v1_01
//
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include "TempSwimSpeedOptim/TempSwimSpeedOptim.hpp"

//probability foraging optimization of energy gain
//	V,U		Vector of speed of predators and prey at each time step
//	K		Food availability for prey, i.e., the obtained reward will be K*(1 + omega*U)
//	C		Vector of metabolic predation cost for predators
//	L		Vector of influence of brightness on the predation rate
//	d		relative density of predator/prey
//	e		relative predation risk of resting prey to foraging ones
//	omega	influence of the speed of prey on the foraging efficiency, i.e., the obtained reward will be K*(1 + omega*U)
//following three parameters determine the predation rate: a*(v-u)^b / {1 + h*a*(v-u)^b}  
//	a		coefficienct of the predation rate (now it is fixed to one)
//	b		non-linear influence of speed difference
//	h		average handling time for predation a prey
//following twp parameters determine the cost of prey
//	cb		metabolic cost for prey (should pay both for resting and foraging)
//	cf		foraging cost for prey (should pay only for foraging)
// [[Rcpp::export]]
Rcpp::List tss_probforage_energygain_optimize(
	Rcpp::NumericVector V,
	Rcpp::NumericVector U,
	Rcpp::NumericVector K,
	Rcpp::NumericVector C,
	Rcpp::NumericVector L,
	double d,
	double e,
	double omega,
	double b,
	double h,
	double cb,
	double cf
){
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::exp_ratio_fitness;
	using predator_fitness = tempss::difference_fitness;
	using this_system = tempss::stepdrdm_optimizer_system<prey_fitness, predator_fitness>;
	using this_state = typename this_system::state;
	using freq_state = tempss::freq_state;

	this_system System(
		predation(1.0, b, h), 
		prey_fitness(1.0, 0.0), 
		predator_fitness(), 
		V.begin(), V.end(), 
		U.begin(), U.end(), 
		K.begin(), K.end(), 
		C.begin(), C.end(), 
		L.begin(), L.end(),
		d, e,omega,cb,cf
	);

	this_state svLower;
	this_state svUpper;
	std::tie(svLower, svUpper) = System();

	freq_state vLower = System.get_freq_state(svLower);
	freq_state vUpper = System.get_freq_state(svUpper);

	double PreyWL = System.get_prey_fitness(svLower);
	double PreyWH = System.get_prey_fitness(svUpper);
	double PredatorWL = System.get_predator_fitness(svLower);
	double PredatorWH = System.get_predator_fitness(svUpper);

	double PreyW = 0.0;
	double PredatorW = 0.0;
	Rcpp::NumericVector Prey(vLower.begin(), vLower.end());
	Rcpp::NumericVector Predator(vLower.size());

	Rcpp::NumericVector Thr(vLower.size());
	Rcpp::NumericVector PreyReward(vLower.size());
	Rcpp::NumericVector PreyCost(vLower.size());
	Rcpp::NumericVector PredatorReward(vLower.size());
	Rcpp::NumericVector PredatorCost(vLower.size());

	// Lower
	if(PreyWL>PreyWH) {
		PreyW = PreyWL;
		PredatorW = PredatorWL;
		for(unsigned int i = 0; i < System.size(); ++i) {
			const auto& TimeInfo = System.at(i);
			Prey[i] = TimeInfo.prey_strategy(svLower[i]);
			Predator[i] = TimeInfo.predator_strategy(svLower[i]);
			Thr[i] = TimeInfo.f_threshold();
			PreyReward[i] = TimeInfo.prey_reward(svLower[i]);
			PreyCost[i] = TimeInfo.prey_cost(svLower[i]);
			PredatorReward[i] = TimeInfo.predator_reward(svLower[i]);
			PredatorCost[i] = TimeInfo.predator_cost(svLower[i]);
		}
	}// Upper
	else {
		PreyW = PreyWH;
		PredatorW = PredatorWH;
		for(unsigned int i = 0; i < System.size(); ++i) {
			const auto& TimeInfo = System.at(i);
			Prey[i] = TimeInfo.prey_strategy(svUpper[i]);
			Predator[i] = TimeInfo.predator_strategy(svUpper[i]);
			Thr[i] = TimeInfo.f_threshold();
			PreyReward[i] = TimeInfo.prey_reward(svUpper[i]);
			PreyCost[i] = TimeInfo.prey_cost(svUpper[i]);
			PredatorReward[i] = TimeInfo.predator_reward(svUpper[i]);
			PredatorCost[i] = TimeInfo.predator_cost(svUpper[i]);
		}
	}

	return Rcpp::List::create(
		Rcpp::Named("Prey") = Prey,
		Rcpp::Named("Predator") = Predator,
		Rcpp::Named("PreyW") = PreyW,
		Rcpp::Named("PredatorW") = PredatorW,
		Rcpp::Named("ThresholdPreyFreq") = Thr,
		Rcpp::Named("PreyReward") = PreyReward,
		Rcpp::Named("PreyCost") = PreyCost,
		Rcpp::Named("PredatorReward") = PredatorReward,
		Rcpp::Named("PredatorCost") = PredatorCost
	);
}

// [[Rcpp::export]]
Rcpp::List tss_probforage_energygain_optimize_linear(
	Rcpp::NumericVector V,
	Rcpp::NumericVector U,
	Rcpp::NumericVector K,
	Rcpp::NumericVector C,
	Rcpp::NumericVector L,
	double d,
	double e,
	double omega,
	double b,
	double h,
	double cb,
	double cf
) {
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::linear_ratio_fitness;
	using predator_fitness = tempss::difference_fitness;
	using this_system = tempss::stepdrdm_optimizer_system<prey_fitness, predator_fitness>;
	using this_state = typename this_system::state;
	using freq_state = tempss::freq_state;

	this_system System(
		predation(1.0, b, h),
		prey_fitness(0.0),
		predator_fitness(),
		V.begin(), V.end(),
		U.begin(), U.end(),
		K.begin(), K.end(),
		C.begin(), C.end(),
		L.begin(), L.end(),
		d, e, omega, cb, cf
	);

	this_state svLower;
	this_state svUpper;
	std::tie(svLower, svUpper) = System();

	freq_state vLower = System.get_freq_state(svLower);
	freq_state vUpper = System.get_freq_state(svUpper);

	double PreyWL = System.get_prey_fitness(svLower);
	double PreyWH = System.get_prey_fitness(svUpper);
	double PredatorWL = System.get_predator_fitness(svLower);
	double PredatorWH = System.get_predator_fitness(svUpper);

	double PreyW = 0.0;
	double PredatorW = 0.0;
	Rcpp::NumericVector Prey(vLower.begin(), vLower.end());
	Rcpp::NumericVector Predator(vLower.size());

	Rcpp::NumericVector Thr(vLower.size());
	Rcpp::NumericVector PreyReward(vLower.size());
	Rcpp::NumericVector PreyCost(vLower.size());
	Rcpp::NumericVector PredatorReward(vLower.size());
	Rcpp::NumericVector PredatorCost(vLower.size());

	// Lower
	if(PreyWL>PreyWH) {
		PreyW = PreyWL;
		PredatorW = PredatorWL;
		for(unsigned int i = 0; i < System.size(); ++i) {
			const auto& TimeInfo = System.at(i);
			Prey[i] = TimeInfo.prey_strategy(svLower[i]);
			Predator[i] = TimeInfo.predator_strategy(svLower[i]);
			Thr[i] = TimeInfo.f_threshold();
			PreyReward[i] = TimeInfo.prey_reward(svLower[i]);
			PreyCost[i] = TimeInfo.prey_cost(svLower[i]);
			PredatorReward[i] = TimeInfo.predator_reward(svLower[i]);
			PredatorCost[i] = TimeInfo.predator_cost(svLower[i]);
		}
	}// Upper
	else {
		PreyW = PreyWH;
		PredatorW = PredatorWH;
		for(unsigned int i = 0; i < System.size(); ++i) {
			const auto& TimeInfo = System.at(i);
			Prey[i] = TimeInfo.prey_strategy(svUpper[i]);
			Predator[i] = TimeInfo.predator_strategy(svUpper[i]);
			Thr[i] = TimeInfo.f_threshold();
			PreyReward[i] = TimeInfo.prey_reward(svUpper[i]);
			PreyCost[i] = TimeInfo.prey_cost(svUpper[i]);
			PredatorReward[i] = TimeInfo.predator_reward(svUpper[i]);
			PredatorCost[i] = TimeInfo.predator_cost(svUpper[i]);
		}
	}

	return Rcpp::List::create(
		Rcpp::Named("Prey") = Prey,
		Rcpp::Named("Predator") = Predator,
		Rcpp::Named("PreyW") = PreyW,
		Rcpp::Named("PredatorW") = PredatorW,
		Rcpp::Named("ThresholdPreyFreq") = Thr,
		Rcpp::Named("PreyReward") = PreyReward,
		Rcpp::Named("PreyCost") = PreyCost,
		Rcpp::Named("PredatorReward") = PredatorReward,
		Rcpp::Named("PredatorCost") = PredatorCost
	);
}
