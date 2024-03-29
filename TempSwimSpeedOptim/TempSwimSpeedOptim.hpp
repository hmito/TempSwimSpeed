#ifndef TEMPSWIMSPEEDOPTIM_INC
#define TEMPSWIMSPEEDOPTIM_INC 100
#
#include<tuple>
#include<vector>
#include<algorithm>
#include<numeric>
#include<boost/math/tools/roots.hpp>
#include<boost/math/tools/minima.hpp>

namespace tempss{
	namespace detail{
		template<typename T, typename eT, typename ans_type = decltype(std::declval<eT>().predation_threshold(1.0, 1.0, 1.0, 1.0))>
		double find_predation_threshold_impl(const T& val, double v, double u, double l, double c, const eT& eval){
			return val.predation_threshold(v, u, l, c);
		}
		template<typename T>
		double find_predation_threshold_impl(const T& val, double v, double u, double l, double c, ...){
			//val should be a monotonically increasing function for f.
			auto Func = [=](double f) {return c - val(f, v, u, l); };
			double v1 = Func(1);
			double v0 = Func(0);
			if(v1 < 0 && v0 < 0)return 0;
			if(v1 > 0 && v0 > 0)return 1;
			auto Ans = boost::math::tools::bisect(Func, 0., 1., boost::math::tools::eps_tolerance<double>(10));
			return (Ans.first + Ans.second) / 2.0;
		}
	}
	template<typename predation_rate>
	double find_predation_threshold(predation_rate&& Predation_, double v, double u, double l, double c){ return detail::find_predation_threshold_impl(Predation_, v, u, l, c, Predation_); }

	namespace detail{
		template<typename T, typename eT, typename ans_type = decltype(std::declval<eT>().is_increasing(1.0, 1.0, 1.0))>
		bool is_prey_fitness_increasing_impl(const T& val, double r, double m, double drdm, const eT& eval){
			return val.is_increasing(r, m, drdm);
		}
		template<typename T>
		bool is_prey_fitness_increasing_impl(const T& val, double r, double m, double drdm, ...){
			constexpr double d = 1e-2;
			return val(r + drdm * d, m + d) > val(r - drdm*d, m - d);
		}
	}
	template<typename prey_fitness>
	bool is_prey_fitness_increasing(prey_fitness PreyFitness, double r, double m, double drdm){ return detail::is_prey_fitness_increasing_impl(PreyFitness, r, m, drdm, PreyFitness); }

	using freq_state = std::vector<double>;

	template<typename prey_fitness_, typename predator_fitness_>
	struct stepdrdm_optimizer_system {
		using prey_fitness = prey_fitness_;
		using predator_fitness = predator_fitness_;
		using state_element = unsigned char;
		using state = std::vector<state_element>;
	private:
		struct time_info{
		private://given parameters
			double v;			//speed of predator
			double u;			//speed of prey
			double k;			//foraging performance of prey, i.e., f*k*u is the obtained reward
			double c;			//metaboric cost of predator
			double l;			//predation rate
			double d;			//relative density of predator/prey
			double e;			//relative risk of predation for resting prey
			double cb;			//metabolic cost for prey (should pay both for resting and foraging)
			double cf;			//foraging cost for prey (should pay only for foraging)
		private://calculated parameters
			double f_thr;		//thrsholf f of predator foraging
			double p0;
			double pF;
			double p1;
			double r0;			//reward under f=0
			double rF;			//reward under f=f_thr
			double r1;			//reward under f=1
			double m0;			//mortality under f=0
			double mF;			//mortality under f=f_thr
			double m1;			//mortality under f=1
		public:
			template<typename predation>
			time_info(const predation& Predation_, double v_, double u_, double k_, double c_, double l_, double d_, double e_, double omega_, double cb_, double cf_)
				: v(v_)
				, u(u_)
				, k(k_)
				, c(c_)
				, l(l_)
				, d(d_)
				, e(e_)
				, cb(cb_)
				, cf(cf_){

				//necessary effective fraction of prey (i.e., f + e*(1-f))
				//F: fraction of attackable prey for predator
				//f: fraction of foraging prey
				//	F = f + e*(1-f)
				//	f = (F-e)/(1-e) 
				double F_thr = find_predation_threshold(Predation_, v, u, l, c);
				f_thr = (F_thr - e) / (1 - e + std::numeric_limits<double>::min());

				if(f_thr <= 0){
					//predator always forage regadless of prey behaviour
					p0 = Predation_(e, v, u, l);
					p1 = Predation_(1, v, u, l);
					m0 = p0 * d + cf * e + cb;
					m1 = p1 * d + cf * 1 + cb;
					r0 = 0;
					r1 = k*(1+omega_*u);

					// xF == x0
					pF = p0;
					mF = m0;
					rF = r0;
				} else if(f_thr >= 1.0){
					//predator always rest regadless of prey behaviour
					p0 = 0.0;
					p1 = 0.0;
					m0 = 0.0 + cf * e + cb;
					m1 = 0.0 + cf * 1 + cb;
					r0 = 0;
					r1 = k*(1+omega_*u);

					// xF == x1
					pF = p1;
					mF = m1;
					rF = r1;
				} else{
					//predator forage only when f > f_thr (i.e., f == 1)
					p0 = 0.0;
					pF = 0.0;
					p1 = Predation_(1, v, u, l);
					m0 = 0.0 +  cf * e + cb;
					mF = 0.0 + cf * F_thr + cb ;
					m1 = p1*d + cf * 1 + cb;
					r0 = 0;
					rF = k*(1+omega_*u)*f_thr;
					r1 = k*(1+omega_*u);
				}
			}
		public:
			double predator_strategy(state_element s)const{
				if(s == 0)return p0>0 ? 1.0 : 0.0;
				else if(s == 1)return pF>0 ? 1.0 : 0.0;
				else return p1>0 ? 1.0 : 0.0;
			}
			double predator_reward(state_element s)const{
				if(s == 0)return p0;
				else if(s == 1)return pF;
				else return p1;
			}
			double predator_cost(state_element s)const{
				return c*predator_strategy(s);
			}
			double prey_strategy(state_element s)const{
				if(s == 0)return 0.0;
				else if(s == 1)return f_thr;
				else return 1.0;
			}
			double prey_reward(state_element s)const{
				if(s == 0)return r0;
				else if(s == 1)return rF;
				else return r1;
			}
			double prey_cost(state_element s)const{
				if(s == 0)return m0;
				else if(s == 1)return mF;
				else return m1;
			}
			bool is_two_step()const{ return (0.0 < f_thr && f_thr < 1.0); }
			double f_threshold()const{ return f_thr; }
		};
		using container = std::vector<time_info >;
	public:
		using iterator = typename container::iterator;
		using const_iterator = typename container::const_iterator;
	private:
		prey_fitness PreyFitness;
		predator_fitness PredatorFitness;
		container Container;
	public:
		template<typename predation, typename iterator>
		stepdrdm_optimizer_system(predation Predation_, prey_fitness PreyFitness_, predator_fitness PredatorFitness_, 
			iterator VBeg, iterator VEnd, 
			iterator UBeg, iterator UEnd, 
			iterator KBeg, iterator KEnd, 
			iterator CBeg, iterator CEnd,
			iterator LBeg, iterator LEnd, 
			double d_, double e_, double omega_, double cb_, double cf_)
			: PreyFitness(std::move(PreyFitness_))
			, PredatorFitness(std::move(PredatorFitness_)){
			if(std::distance(VBeg, VEnd) != std::distance(UBeg, UEnd)
				|| std::distance(VBeg, VEnd) != std::distance(KBeg, KEnd)
				|| std::distance(VBeg, VEnd) != std::distance(CBeg, CEnd)
				|| std::distance(VBeg, VEnd) != std::distance(LBeg, LEnd)
				){
				throw std::exception();
			}

			while(VBeg != VEnd) {
				Container.emplace_back(Predation_, *VBeg++, *UBeg++, *KBeg++, *CBeg++, *LBeg++, d_, e_, omega_, cb_, cf_);
			}
		}
		std::pair<state,state> operator()(void)const{
			std::vector<double> Data;

			//Add data
			for(const auto& Info : Container){
				if(Info.is_two_step()){
					//predator's forage only when prey fully use this time
					//	set difference between no-use and partial-use 
					Data.push_back(
						(Info.prey_reward(1) - Info.prey_reward(0)) / std::max(Info.prey_cost(1) - Info.prey_cost(0), std::numeric_limits<double>::min())
					);
				} else{
					//predator's behaviour is independent from prey's behaviour
					//	set difference between no-use and full-use
					Data.push_back(
						(Info.prey_reward(2) - Info.prey_reward(0)) / std::max(Info.prey_cost(2) - Info.prey_cost(0), std::numeric_limits<double>::min())
					);
				}
			}

			state Lower(Container.size(), 0);
			state Upper(Container.size(), 0);

			auto No = std::distance(Data.begin(), std::max_element(Data.begin(), Data.end()));
			while(Data[No]>0.0){
				auto& Info = Container.at(No);

				//Update Strategy
				if(Info.is_two_step() && Upper[No]==0){
					//change from no-use to partial-use
					Upper[No] = 1;
					//set difference between partial-use and full-use 
					Data[No] = (Info.prey_reward(2) - Info.prey_reward(1)) / (Info.prey_cost(2) - Info.prey_cost(1) + std::numeric_limits<double>::min());
				} else{
					//change from no-use/partial-use to full-use
					Upper[No] = 2;
					//set drdm to zero (i.e., no reward)
					Data[No] = 0.0;
				}

				//Find next best step No
				auto NewNo = std::distance(Data.begin(), std::max_element(Data.begin(), Data.end()));

				//Check if Fitness is still increasing
				double TotalR = 0.0;
				double TotalM = 0.0;
				for(unsigned int i = 0; i < Container.size(); ++i){
					TotalR += Container.at(i).prey_reward(Upper[i]);
					TotalM += Container.at(i).prey_cost(Upper[i]);
				}

				if(!is_prey_fitness_increasing(PreyFitness, TotalR, TotalM, Data[NewNo])){
					break;
				}

				Lower[No] = Upper[No];
				No = NewNo;
			}

			return std::make_pair(std::move(Lower), std::move(Upper));
		}
		const_iterator begin()const { return std::begin(Container); }
		const_iterator end()const { return std::end(Container); }
		const time_info& at(unsigned int i)const{ return Container.at(i); }
		unsigned int size()const{ return Container.size(); }
		double get_prey_fitness(const state& x){
			double R = 0;
			double M = 0;
			for (unsigned int i = 0; i < x.size(); ++i) {
				R += Container[i].prey_reward(x[i]);
				M += Container[i].prey_cost(x[i]);
			}

			return PreyFitness(R, M);
		}
		double get_predator_fitness(const state& x){
			double R = 0;
			double M = 0;
			for(unsigned int i = 0; i < x.size(); ++i){
				R += Container[i].predator_reward(x[i]);
				M += Container[i].predator_cost(x[i]);
			}

			return PredatorFitness(R, M);
		}
		freq_state get_freq_state(const state& x){
			freq_state y;
			for(unsigned int i = 0; i < Container.size(); ++i){
				y.push_back(Container.at(i).prey_strategy(x[i]));
			}
			return y;
		}
	};

	//Perdation Assumption C1: alpha*(v-u)^beta / {1 + a*(v-u)^b*h} 
	struct probforage_predation{
	private:
		double a;	//inverse of searching time
		double b;	//non-linear influence of speed difference
		double h;	//average handling time for predation a prey
	public:
		probforage_predation(double a_, double b_, double h_)
			: a(a_)
			, b(b_)
			, h(h_){
		}
		double operator()(double pf, double v, double u, double l) const{
			//when beta = 0, return value is non-zero constant when v < u.
			//when beta > 0, return value is zero when v<u.
			double gamma =  a*l*pf*(b > 0 ? std::pow(std::max(v - u,0.0), b):1.0);
			return gamma / (1 + gamma*h);
		}
	};

	//Perdation Assumption C1_fix: alpha*(v-u)^beta / {1 + a*(v-u)^b*h} 
	struct probforage_fix_predation{
	private:
		double a;	//inverse of searching time
		double b;	//non-linear influence of speed difference
		double h;	//average handling time for predation a prey
	public:
		probforage_fix_predation(double a_, double b_, double h_)
			: a(a_)
			, b(b_)
			, h(h_){
		}
		double operator()(double pf, double v, double u, double l) const{
			//return value is always zero when v<u.
			double gamma =  0.0;
			if(v>u) gamma = a*l*pf*std::pow(v - u, b);
			return gamma / (1 + gamma*h);
		}
	};

	//Perdation Assumption C4: maxr / [1+exp(-beta*(v-u-dv50))]
	struct sigmoid_predation{
	private:
		double maxr;	//non-linear influence of speed difference
		double beta;	//average handling time for predation a prey
		double dvh; 	//inverse of searching time
	public:
		sigmoid_predation(double maxr_, double beta_, double dvh_)
			: maxr(maxr_)
			, beta(beta_)
			, dvh(dvh_){
		}
		double operator()(double pf, double v, double u, double l) const{
			return l * pf * maxr / (1+ std::exp(-beta*(v-u-dvh)));
		}
	};

	//Perdation Assumption C2: 1 / {(s/v/f + d/(v-u) + h} 
	struct search_hunt_forage_predation {
	private:
		double s;	//average search time when the prey density is 1 (full predation) and swim speed is 1
		double d;	//average distance that the prey start to run away
		double h;	//average handling time for predation a prey
	public:
		search_hunt_forage_predation(double s_, double d_, double h_)
			: s(s_)
			, d(d_)
			, h(h_) {
		}
		double operator()(double f, double v, double u) const {
			if (f <= 0 || v <= 0 || v <= u )return 0.;
			return 1 / (s / (v*f) + d / (v - u) + h);
		}
	};

	//Perdation Assumption C3: 1 / {(s/v/f/exp(omega*(v-u)) + h)}
	struct search_probforage_predation {
	private:
		double omega;	//function of successful hunting for given relative speed
		double s;		//average search time when the prey density is 1 (full predation) and swim speed is 1
		double h;		//average handling time for predation a prey
	public:
		search_probforage_predation(double omega_, double s_, double h_)
			: omega(omega_)
			, s(s_)
			, h(h_) {
		}
		double operator()(double f, double v, double u) const {
			if (f <= 0 || v <= 0)return 0.;
			if(v <= u)return 0.;
			double lambda = 1 - std::exp(omega*(v - u));
			return 1 / (s / v / f / lambda + h);
		}
	};

	//fitness reward-mortality ratio with basic mortality rate: {1 - alpha*exp(-beta*r)} / {m + base_m}
	struct exp_ratio_fitness{
	private:
		double alpha;
		double base_m;
	public:
		exp_ratio_fitness(double alpha_, double base_m_) :alpha(alpha_), base_m(base_m_){}
		double operator()(double r, double m) const{
			return (1 - std::exp(-alpha*r)) / (m + base_m);
		}
		bool is_increasing(double r, double m, double drdm)const{
			if(drdm > 1.0e100)return true;
			return (base_m+m)*alpha*std::exp(-alpha*r)*drdm > 1 - std::exp(-alpha*r);
		}
		double minimum_drdm(double r, double m)const{
			return (1 - std::exp(-alpha*r)) / ((base_m + m)*alpha*std::exp(-alpha*r) + std::numeric_limits<double>::min());
		}
	};

	//fitness reward-mortality ratio with basic mortality rate: {1 - alpha*exp(-beta*r)} / {m + base_m}
	struct linear_ratio_fitness{
	private:
		double base_m;
	public:
		linear_ratio_fitness(double base_m_) :base_m(base_m_){}
		double operator()(double r, double m) const{
			return r / (base_m + m);
		}
		bool is_increasing(double r, double m, double drdm)const{
			return (base_m + m)*drdm > r;
		}
		double minimum_drdm(double r, double m)const{
			return r / (base_m + m);
		}
	};

	//fitness reward-mortality ratio with basic mortality rate: {1 - alpha*exp(-beta*r)} / {m + base_m}
	struct difference_fitness{
	public:
		double operator()(double r, double m) const{
			return r - m;
		}
		bool is_increasing(double r, double m, double drdm)const{
			return drdm > 1;
		}
		double minimum_drdm(double r, double m)const{
			return 1.0;
		}
	};
}
#
#endif
