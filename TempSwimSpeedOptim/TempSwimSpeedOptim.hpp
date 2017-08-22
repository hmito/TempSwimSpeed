#ifndef TEMPSWIMSPEEDOPTIM_INC
#define TEMPSWIMSPEEDOPTIM_INC 100
#
#include<tuple>
#include<vector>
#include<algorithm>
#include<numeric>
#include<boost/math/tools/roots.hpp>
#include<boost/math/tools/minima.hpp>
#include"../../hmLib/optimize.hpp"
#include"../../hmLib/random.hpp"

namespace hmLib {
	namespace optimize {
		/*!
		@brief	golden section search method
		Search the maximum value of Func in [MinVal:MaxVal].
		@param [in]	Func	Target function. The function need not to have a range where the slope becomes zero.
		@param [in]	MinVal	Minimum value of searching range.
		@param [in]	MaxVal	Maximum value of searching range.
		@param [in]	Error	Accuracy of the search. The return range will be less than Error.
		@return	Pair of the value (Min, Max) where Func becomes maximum value.
				The distance of return pair Max-Min will be less than Error if the searching is sucessfully finished.
				If the searching of maximum value is failed, the return pair distance will be larger than Error.*/
		template<typename func, typename state>
		std::pair<state, state> golden_section_search(func&& Func, state MinVal, state MaxVal, double Error) {
			constexpr double gratio = 1.6180339887498948482045868343656;

			auto MinEval = Func(MinVal);
			auto MaxEval = Func(MaxVal);

			//First division: MinSide
			auto LowerVal = (gratio*MinVal + MaxVal) / (1 + gratio);
			auto LowerEval = Func(LowerVal);
			if (LowerEval < MinEval && LowerEval < MaxEval)return std::make_pair(MinVal, MaxVal);

			//Second division: MaxSide
			state UpperVal = (gratio*LowerVal + MaxVal) / (1 + gratio);
			auto UpperEval = Func(UpperVal);

			while (MaxVal - MinVal > Error) {
				//Update Min 
				if (LowerEval < UpperEval) {
					MinVal = std::move(LowerVal);
					MinEval = std::move(LowerEval);

					LowerVal = std::move(UpperVal);
					LowerEval = std::move(UpperEval);

					UpperVal = (gratio*LowerVal + MaxVal) / (1 + gratio);
					UpperEval = Func(UpperVal);
				}
				//Update Max
				else{
					MaxVal = std::move(UpperVal);
					MaxEval = std::move(UpperEval);

					UpperVal = std::move(LowerVal);
					UpperEval = std::move(LowerEval);

					LowerVal = (MinVal + gratio*UpperVal) / (1 + gratio);
					LowerEval = Func(LowerVal);
				}
			}

			return std::make_pair(MinVal, MaxVal);
		}

		/*!
		@brief	golden section search method with small division method.
		Search the maximum value of Func in [MinVal:MaxVal].
		@param [in]	Func	Target function. The function can have a range where the slope becomes zero.
		@param [in]	MinVal	Minimum value of searching range.
		@param [in]	MaxVal	Maximum value of searching range.
		@param [in]	Error	Accuracy of the search. The return range will be less than Error.
		@return	Pair of the value (Min, Max) where Func becomes maximum value.
				The distance of return pair Max-Min will be less than Error if the searching is sucessfully finished.
				If the searching of maximum value is failed, the return pair distance will be larger than Error.*/
		template<typename func, typename state>
		std::pair<state, state> flatable_golden_section_search(func&& Func, state MinVal, state MaxVal, double Error) {
			constexpr double gratio = 1.6180339887498948482045868343656;

			auto MinEval = Func(MinVal);
			auto MaxEval = Func(MaxVal);

			//First division: MinSide
			auto LowerVal = (gratio*MinVal + MaxVal) / (1 + gratio);
			auto LowerEval = Func(LowerVal);
			if (LowerEval < MinEval && LowerEval < MaxEval)return std::make_pair(MinVal, MaxVal);

			//Second division: MaxSide
			state UpperVal = (gratio*LowerVal + MaxVal) / (1 + gratio);
			auto UpperEval = Func(UpperVal);

			while (MaxVal - MinVal > Error) {
				//Update Min 
				if (LowerEval < UpperEval) {
					MinVal = std::move(LowerVal);
					MinEval = std::move(LowerEval);

					LowerVal = std::move(UpperVal);
					LowerEval = std::move(UpperEval);

					UpperVal = (gratio*LowerVal + MaxVal) / (1 + gratio);
					UpperEval = Func(UpperVal);
				}
				//Update Max
				else if (UpperEval < LowerEval) {
					MaxVal = std::move(UpperVal);
					MaxEval = std::move(UpperEval);

					UpperVal = std::move(LowerVal);
					UpperEval = std::move(LowerEval);

					LowerVal = (MinVal + gratio*UpperVal) / (1 + gratio);
					LowerEval = Func(LowerVal);
				}
				//Same Case
				else {
					unsigned int Div = 11;

					unsigned int BestNoLower;
					unsigned int BestNoUpper;

					while (true) {
						BestNoLower = 0;
						BestNoUpper = 0;
						auto BestEval = MinEval;
						for (unsigned int i = 1; i < Div - 1; ++i) {
							state Val = (MinVal*(Div - i) + MaxVal*i) / Div;
							auto Eval = Func(Val);

							if (BestEval < Eval) {
								BestNoLower = i;
								BestNoUpper = i;
								BestEval = Eval;
							} else if (BestEval == Eval) {
								BestNoUpper = i;
							}
						}
						if (BestEval < MaxEval) {
							BestNoLower = Div;
							BestNoUpper = Div;
							BestEval = MaxEval;
						} else if (BestEval == MaxEval) {
							BestNoUpper = Div;
						}

						if (BestNoLower != 0 || BestNoUpper != Div)break;

						if ((MaxVal - MinVal) / Div < Error) {
							return std::make_pair(MinVal, MaxVal);
						}
						Div = (Div - 1) * 10 + 1;
					}
					if (BestNoLower > 0)BestNoLower -= 1;
					if (BestNoUpper < Div)BestNoUpper += 1;
					auto NewMinVal = (MinVal*(Div - BestNoLower) + MaxVal*BestNoLower) / Div;
					auto NewMaxVal = (MinVal*(Div - BestNoUpper) + MaxVal*BestNoUpper) / Div;

					if (NewMinVal <= MinVal && NewMaxVal >= MaxVal) {
						return std::make_pair(MinVal, MaxVal);
					}
					MinVal = NewMinVal;
					MaxVal = NewMaxVal;
					MinEval = Func(MinVal);
					MaxEval = Func(MaxVal);

					//First division: MinSide
					auto LowerVal = (gratio*MinVal + MaxVal) / (1 + gratio);
					auto LowerEval = Func(LowerVal);
					if (LowerEval < MinEval && LowerEval < MaxEval)return std::make_pair(MinVal, MaxVal);

					//Second division: MaxSide
					state UpperVal = (gratio*LowerVal + MaxVal) / (1 + gratio);
					auto UpperEval = Func(UpperVal);
				}
			}

			return std::make_pair(MinVal, MaxVal);
		}
	}
}
namespace tempss{
	namespace detail {
		template<typename T, typename eT, typename ans_type = decltype(std::declval<eT>().predation_threshold(1.0, 1.0, 1.0))>
		double find_predation_threshold_impl(const T& val, double v, double u, double c, const eT& eval) {
			return val.predation_threshold(v, u, c);
		}
		template<typename T>
		double find_predation_threshold_impl(const T& val, double v, double u, double c, ...) {
			auto Func = [=](double f) {return c - f*val(f, v, u); };
			double v1 = Func(1);
			double v0 = Func(0);
			if (v1 < 0 && v0 < 0)return 0;
			if (v1 > 0 && v0 > 0)return 1;
			auto Ans = boost::math::tools::bisect(Func, 0., 1., boost::math::tools::eps_tolerance<double>(10));
			return (Ans.first + Ans.second) / 2.0;
		}
	}
	template<typename predation_rate>
	double find_predation_threshold(predation_rate&& Predation_, double v, double u, double c) { return detail::find_predation_threshold_impl(Predation_, v, u, c, Predation_); }

	namespace detail {
		template<typename T, typename eT, typename ans_type = decltype(std::declval<eT>().is_increasing(1.0, 1.0, 1.0))>
		bool is_prey_fitness_increasing_impl(const T& val, double r, double m, double drdm, const eT& eval) {
			return val.is_increasing(r, m, drdm);
		}
		template<typename T>
		bool is_prey_fitness_increasing_impl(const T& val, double r, double m, double drdm, ...) {
			constexpr double d = 1e-2;
			return val(r + drdm * d, m + d) > val(r - drdm*d, m - d);
		}
	}
	template<typename prey_fitness>
	bool is_prey_fitness_increasing(prey_fitness PreyFitness, double r, double m, double drdm) { return detail::is_prey_fitness_increasing_impl(PreyFitness, r, m, drdm, PreyFitness); }

	template<typename predation>
	struct time_info {
	private:
		//any function of
		//	predation rate of predator (predation probability of prey)
		const predation& Predation;
	private://given parameters
		double v;			//speed of predator
		double u;			//speed of prey
		double mu;			//basic mortality change of prey by foraging
		double d;			//relative density of predator/prey
	private://calculated parameters
		double c;			//metaboric cost of predator
		double r;			//predation perfodrdmance of prey
		double f_thr;		//thrsholf of predator foragingif(
		double drdm0;		//drdm at f=0
		double drdm1;		//drdm at f=1
		double drdmF;		//drdm at f=f_thr+delta
	public:
		template<typename prey_reward_t, typename predator_cost_t>
		time_info(const predation& Predation_, const prey_reward_t& PreyReward, const predator_cost_t& PredatorCost, double v_, double u_, double mu_, double d_)
			: Predation(Predation_)
			, v(v_)
			, u(u_)
			, mu(mu_)
			, d(d_) {
			c = PredatorCost(v);
			r = PreyReward(u);
			f_thr = find_predation_threshold(Predation_, v, u, c);
			if (f_thr <= 0.0)f_thr = -1e-10;
			drdm0 = std::max(0., r / (prey_mortality(0) + std::numeric_limits<double>::min()));
			drdm1 = std::max(0., r / (prey_mortality(1) + std::numeric_limits<double>::min()));
			drdmF = f_thr <= 0 ? drdm0 : f_thr >= 1 ? drdm1 : std::max(0., r / (prey_mortality(f_thr, 1) + std::numeric_limits<double>::min()));
		}
	public:
		double predator_strategy(double f)const {
			return f > f_thr? 1 : 0;
		}
		double predator_payoff(double f)const {
			return std::max(Predation(f, v, u) - c, 0.0);
		}
		double prey_reward(double f)const {
			return r;
		}
		double prey_mortality(double f)const {
			return prey_mortality(f, predator_strategy(f));
		}
		double prey_mortality(double f, double p) const {
			return mu + Predation(f, v, u) * d * p / (f + std::numeric_limits<double>::min());
		}
		double prey_predation_mortality(double f)const{
			return prey_predation_mortality(f, predator_strategy(f));
		}
		double prey_predation_mortality(double f, double p) const{
			return Predation(f, v, u) * d * p / (f + std::numeric_limits<double>::min());
		}
		double f_threshold()const{ return f_thr; }
		double prey_drdm0()const { return drdm0; }
		double prey_drdm1()const{ return drdm1; }
		double prey_drdmF()const{ return drdmF; }
	};

	using state = std::vector<double>;

	template<typename info_iterator, typename fitness>
	std::pair<state, state> stepdrdm_bisect_optimize(info_iterator InfoBeg, info_iterator InfoEnd, fitness&& Fitness){
		using data_t = std::tuple<unsigned int, double, double>;
		std::vector<data_t> Data;

		//Add data
		for(auto Itr = InfoBeg; Itr != InfoEnd; ++Itr){
			double f = Itr->f_threshold();
			unsigned int Pos = static_cast<unsigned int>(std::distance(InfoBeg, Itr));
			if(f < 1.0){
				Data.emplace_back(Pos, Itr->prey_drdm0(), f);
				Data.emplace_back(Pos, Itr->prey_drdm1(), 1.0);
			} else{
				Data.emplace_back(Pos, Itr->prey_drdm1(), 1.0);
			}
		}

		std::sort(Data.begin(), Data.end(), [](const data_t& v1, const data_t& v2){return std::get<1>(v1) > std::get<1>(v2); });

		unsigned int Size = static_cast<unsigned int>(std::distance(InfoBeg, InfoEnd));
		std::vector<double> Lower(Size, 0.);
		std::vector<double> Upper(Size, 0.);

		unsigned int i = 0;
		for(; i < Data.size(); ++i){
			unsigned int No = std::get<0>(Data[i]);
			double drdm = std::get<1>(Data[i]);
			double f = std::get<2>(Data[i]);

			Upper[No] = f;

			double r = 0.0;
			double m = 0.0;
			for(unsigned int j = 0; j < Size; ++j){
				r += InfoBeg[j].prey_reward(Upper[j])*Upper[j];
				m += InfoBeg[j].prey_mortality(Upper[j])*Upper[j];
			}
			if(!is_prey_fitness_increasing(Fitness, r, m, drdm))break;

			Lower[No] = f;
		}

		return std::make_pair(std::move(Lower), std::move(Upper));
	}

	namespace detail {
		struct mutate {
			double Sigma;
			mutate(double Sigma_):Sigma(Sigma_){}
			template<typename urbg>
			void operator()(state& x, urbg&& RandEng) {
				unsigned int pos = std::uniform_int_distribution<unsigned int>(0, x.size()-1)(RandEng);

				x[pos] = std::max(0.,std::min(1.,std::normal_distribution<double>(x[pos], Sigma)(RandEng)));
			}
		};
	}
	template<typename info_iterator, typename fitness>
	state hill_climbing_search(info_iterator InfoBeg, info_iterator InfoEnd, fitness&& Fitness, unsigned int StableStep, unsigned int MaxStep, double Sigma) {
		state State(std::distance(InfoBeg, InfoEnd), 0.5);
		auto Evaluate = [=,&Fitness](const state& x) {
			double R = 0;
			double M = 0;
			auto Itr = InfoBeg;
			for (unsigned int i = 0; i < x.size(); ++i, ++Itr) {
				R += x[i] * Itr->prey_reward(x[i]);
				M += x[i] * Itr->prey_mortality(x[i]);
			}

			return Fitness(R, M);
		};
		detail::mutate Mutate(Sigma);
		hmLib::optimize::hill_climbing_search(State, Evaluate, Mutate, hmLib::optimize::state_breakers::limited_const_stable_breaker<state, double>(StableStep, MaxStep),hmLib::random::default_engine());
		return State;
	}

	template<typename predation_, typename prey_fitness_>
	struct optimizer_system {
		using predation = predation_;
		using prey_fitness = prey_fitness_;
	private:
		using this_time_info = time_info<predation>;
		using container = std::vector<this_time_info >;
	public:
		using iterator = typename container::iterator;
		using const_iterator = typename container::const_iterator;
	private:
		predation Predation;
		prey_fitness PreyFitness;
		container Container;
	public:
		template<typename prey_reward, typename predator_cost, typename iterator>
		optimizer_system(predation Predation_, prey_fitness PreyFitness_, prey_reward&& PreyReward, predator_cost&& PredatorCost, iterator VBeg, iterator VEnd, iterator UBeg, iterator UEnd, double base_mu_, double d_)
			: Predation(std::move(Predation_))
			, PreyFitness(std::move(PreyFitness_)) {
			if (std::distance(VBeg, VEnd) != std::distance(UBeg, UEnd))throw std::exception();

			for (; VBeg != VEnd; ++VBeg, ++UBeg) {
				Container.emplace_back(Predation, PreyReward, PredatorCost, *VBeg, *UBeg, base_mu_, d_);
			}
		}
		std::pair<state,state> optimize_by_stepbisect(void)const{
			return stepdrdm_bisect_optimize(Container.begin(), Container.end(), PreyFitness);
		}
		state optimize_by_hill_climbing(unsigned int StableStep, unsigned int MaxStep, double Sigma)const{
			return hill_climbing_search(Container.begin(), Container.end(), PreyFitness, StableStep, MaxStep, Sigma);
		}
		const_iterator begin()const { return std::begin(Container); }
		const_iterator end()const { return std::end(Container); }
		const this_time_info& at(unsigned int i)const{ return Container.at(i); }
		unsigned int size()const{ return Container.size(); }
		double get_prey_fitness(const state& x){
			double R = 0;
			double M = 0;
			for (unsigned int i = 0; i < x.size(); ++i) {
				R += x[i] * Container[i].prey_reward(x[i]);
				M += x[i] * Container[i].prey_mortality(x[i]);
			}

			return PreyFitness(R, M);
		}
		double get_predator_fitness(const state& x){
			double W = 0.0;

			for(unsigned int i = 0; i < x.size(); ++i){
				W += Container[i].predator_payoff(x[i]);
			}
			
			return W;
		}
	};

	//Assumption C1: alpha*(v-u)^beta / {1 + a*(v-u)^b*h} 
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
		double operator()(double f, double v, double u) const{
			if(f <= 0 || v <= 0 || v<=u)return 0.;

			double gamma = a*std::pow(v - u, b);
			return gamma * f / (1 + gamma*f*h);
		}
	};

	//Assumption C2: 1 / {(s/v/f + d/(v-u) + h} 
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
			if (f <= 0 || v <= 0)return 0.;
			return 1 / (s / v / f + d / (v - u) + h);
		}
	};

	//Assumption C3: 1 / {(s/v/f/exp(omega*(v-u)) + h)}
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

	//prey reward linearly increase with the prey speed
	struct linear_prey_reward {
	private:
		double k;
	public:
		linear_prey_reward(double k_) :k(k_) {}
		double operator()(double u) const { return k*u; }
	};

	//predator foraging cost is constant
	struct constant_predator_cost {
	private:
		double c;
	public:
		constant_predator_cost(double c_) :c(c_) {}
		double operator()(double v) const { return c; }
	};

	//fitness exponentially saturate
	struct exp_diff_fitness {
	private:
		double lambda;	//speed of saturation
	public:
		exp_diff_fitness(double lambda_) :lambda(lambda_) {}
		double operator()(double r, double m) const {
			return 1 - std::exp(-r) - lambda*m;
		}
	};

	//fitness reward-mortality ratio with basic mortality rate: {1 - alpha*exp(-beta*r)} / {m + base_m}
	struct ratio_fitness{
	private:
		double alpha;
		double base_m;
	public:
		ratio_fitness(double alpha_, double base_m_) :alpha(alpha_), base_m(base_m_){}
		double operator()(double r, double m) const{
			return (1 - std::exp(-alpha*r)) / (m + base_m);
		}
		bool is_increasing(double r, double m, double drdm)const{
			return (base_m+m)*alpha*std::exp(-alpha*r)*drdm > 1 - std::exp(-alpha*r);
		}
	};
}
#
#endif
