#include<vector>
#include<iostream>
#include"TempSwimSpeedOptim.hpp"

int main(void) {
	std::vector<double> V1{ 2.22498176022918,2.1896907799223,2.15154879672934,2.11315512377216,2.07712622641425,2.04591741446998,2.02165551689834,2.00599394192418,2,2.00408216836248,2.01796225398381,2.04069435196683,2.07072930740057,2.10602028770733,2.1441622709003,2.18255594385749,2.21858484121542,2.24979365315971,2.27405555073135,2.28971712570552,2.29571106762972,2.29162889926725,2.27774881364592,2.25501671566291 };
	std::vector<double> U1{ 1.12380720067661,1.04222833053908,1,1,1.04222833053908,1.12380720067661,1.23917714522818,1.38047589313045,1.5380741682211,1.70123190849616,1.85883018358681,2.00012893148908,2.11549887604065,2.19707774617819,2.23930607671726,2.23930607671726,2.19707774617819,2.11549887604065,2.00012893148908,1.85883018358681,1.70123190849616,1.5380741682211,1.38047589313045,1.23917714522818 };
	std::vector<double> V2{ 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 };
	std::vector<double> U2{ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	std::vector<double> K{ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	std::vector<double> C{ 0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02 };
	std::vector<double> L{ 0.991481352436089,0.92670456206338,0.813306973345948,0.676218810601917,0.539389913660382,0.419172056302892,0.322863323548505,0.25090436400304,0.200135342458379,0.16640123307912,0.146039297503055,0.136498062120978,0.136498062120978,0.146039297503055,0.16640123307912,0.200135342458379,0.25090436400304,0.322863323548505,0.419172056302892,0.539389913660382,0.676218810601917,0.813306973345948,0.92670456206338,0.991481352436089 };
	double d = 0.1;
	double e = 0.15;
	double omega = 0.0;
	double b = 0.0;
	double h = 2;
	double cb = 0.01;
	double cf = 0.05;
	using predation = tempss::probforage_predation;
	using prey_fitness = tempss::linear_ratio_fitness;
	using predator_fitness = tempss::difference_fitness;
	using this_system = tempss::stepdrdm_optimizer_system<prey_fitness, predator_fitness>;
	using this_state = typename this_system::state;
	using freq_state = tempss::freq_state;

	this_system System1(
		predation(1.0, b, h),
		prey_fitness(0.0),
		predator_fitness(),
		V1.begin(), V1.end(),
		U1.begin(), U1.end(),
		K.begin(), K.end(),
		C.begin(), C.end(),
		L.begin(), L.end(),
		d, e, omega, cb, cf
	);

	this_system System2(
		predation(1.0, b, h),
		prey_fitness(0.0),
		predator_fitness(),
		V2.begin(), V2.end(),
		U2.begin(), U2.end(),
		K.begin(), K.end(),
		C.begin(), C.end(),
		L.begin(), L.end(),
		d, e, omega, cb, cf
	);

	this_state svLower1;
	this_state svUpper1;
	std::tie(svLower1, svUpper1) = System1();

	this_state svLower2;
	this_state svUpper2;
	std::tie(svLower2, svUpper2) = System2();

	system("pause");
	return 0;
}
