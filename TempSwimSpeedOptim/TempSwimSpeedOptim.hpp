#ifndef TEMPSWIMSPEEDOPTIM_INC
#define TEMPSWIMSPEEDOPTIM_INC
#
#include<tuple>
#include<vector>
#include<algorithm>
#include<numeric>
#include<boost/math/tools/roots.hpp>
#include<boost/math/tools/minima.hpp>
#include"../../hmLib/hmLib/optimize.hpp"
#include"../../hmLib/hmLib/random.hpp"

namespace hmLib {
	namespace optimize {
		//簡単な下に凸な関数f(x) = (x-1)^2。x=1で最小値0。
		double f(double x) {
			return (x - 1.0)*(x - 1.0);
		}

		//黄金分割法
		// f(x)が区間[lb,ub]で凸ならば、その極値を返す
		// 反復時に値が使いまわせるので、fの計算が1回のみでよい
		// ub: 下限    ub: 上限    K: 反復回数
		template<typename func, typename state>
		std::pair<state, state> golden_section_search(func&& Func, state MinVal, state MaxVal, unsigned int Bits) {
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

			for (unsigned int i = 0; i < Bits; ++i) {
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
				else {
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
	}
}
namespace tempss{
	//黄金分割法
	// f(x)が区間[lb,ub]で凸ならば、その極値を返す
	// 反復時に値が使いまわせるので、fの計算が1回のみでよい
	// ub: 下限    ub: 上限    K: 反復回数
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
			else if (UpperEval < LowerEval) {
				MaxVal = std::move(UpperVal);
				MaxEval = std::move(UpperEval);

				UpperVal = std::move(LowerVal);
				UpperEval = std::move(LowerEval);

				LowerVal = (MinVal + gratio*UpperVal) / (1 + gratio);
				LowerEval = Func(LowerVal);
			}
			//Same case
			else {
				LowerVal = (MinVal + LowerVal) / 2.;
				LowerEval = Func(LowerVal);
				UpperVal = (MaxVal + UpperVal) / 2.;
				UpperEval = Func(UpperVal);

				if (LowerVal - MinVal < Error && MaxVal - UpperVal < Error)break;
			}
		}

		return std::make_pair(MinVal, MaxVal);
	}

	namespace detail {
		template<typename T, typename eT, typename ans_type = decltype(std::declval<eT>().predation_threshold(1.0, 1.0, 1.0))>
		double find_predation_threshold_impl(const T& val, double v, double u, double c, const eT& eval) {
			return val.predation_threshold(v, u, c);
		}
		template<typename T>
		double find_predation_threshold_impl(const T& val, double v, double u, double c, ...) {
			auto Func = [=](double f) {return c - val(f, v, u); };
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
		double base_mu;		//basic mortality for prey
		double d;			//relative density of predator/prey
	private://calculated parameters
		double c;			//metaboric cost of predator
		double r;			//predation perfodrdmance of prey
		double Max_drdm;	//max dr/dm
		double Min_drdm;	//min dr/dm
		double f_thr;		//thrsholf of predator foragingif(
	public:
		template<typename prey_reward_t, typename predator_cost_t>
		time_info(const predation& Predation_, const prey_reward_t& PreyReward, const predator_cost_t& PredatorCost, double v_, double u_, double base_mu_, double d_)
			: Predation(Predation_)
			, v(v_)
			, u(u_)
			, base_mu(base_mu_)
			, d(d_) {
			c = PredatorCost(v);
			r = PreyReward(u);
			f_thr = find_predation_threshold(Predation_, v, u, c);
			if (f_thr <= 0.0)f_thr = -1e-10;
			Max_drdm = std::max(0., prey_reward(0) / (prey_mortality(0) + std::numeric_limits<double>::min()));
			Min_drdm = std::max(0., prey_reward(1) / (prey_mortality(1) + std::numeric_limits<double>::min()));
		}
	public:
		double predator_strategy(double f)const {
			return f > f_thr ? 1 : 0;
		}
		double predator_payoff(double f)const {
			return std::max(Predation(f, v, u) - c, 0);
		}
		double prey_reward(double f)const {
			return r;
		}
		double prey_mortality(double f)const {
			return prey_mortality(f, predator_strategy(f));
		}
		double prey_mortality(double f, double p) const {
			return base_mu + Predation(f, v, u) * d * p / (f + std::numeric_limits<double>::min());
		}
		double find_f(double drdm, unsigned int Bits)const {
			if (drdm >= max_drdm()) return 0.0;
			else if (drdm <= min_drdm()) return 1.0;
			auto val = boost::math::tools::bisect(
				[this, drdm](double f) {return prey_reward(f) / (prey_mortality(f) + std::numeric_limits<double>::min()) - drdm; },
				0.0,1.0,
				boost::math::tools::eps_tolerance<double>(Bits)
			);
			return (val.first + val.second) / 2.0;
		}
		double max_drdm()const { return Max_drdm; }
		double min_drdm()const { return Min_drdm; }
		double f_threshold()const { return f_thr; }
	};

	using state = std::vector<double>;

	template<typename info_iterator,typename fitness>
	std::pair<state, state> drdm_bisect_optimize(info_iterator InfoBeg, info_iterator InfoEnd, fitness&& Fitness) {
		double MinRM = std::numeric_limits<double>::max();
		double MaxRM = std::numeric_limits<double>::lowest();
		for (auto Itr = InfoBeg; Itr != InfoEnd; ++Itr) {
			MaxRM = std::max(MaxRM, Itr->max_drdm());
			MinRM = std::min(MinRM, Itr->min_drdm());
		}

		auto Func = [=, &PreyFitness=Fitness](double drdm) {
			double r = 0;
			double m = 0;
			for (auto Itr = InfoBeg; Itr != InfoEnd; ++Itr) {
				double f = Itr->find_f(drdm, 20);
				double dr = Itr->prey_reward(f);
				double dm = Itr->prey_mortality(f);

				r += dr*f;
				m += dm*f;
			}
			return PreyFitness(r, m);
		};

		//Pair first: threshold value, second: 
		auto drdmPair = golden_section_search(Func, MinRM, MaxRM, 1e-10);

		unsigned int Size = std::distance(InfoBeg, InfoEnd);
		std::vector<double> Lower(Size, 0.);
		std::vector<double> Upper(Size, 0.);

		auto Itr = InfoBeg;
		for (unsigned int i = 0; i<Size; ++i, ++Itr) {
			Upper.at(i) = Itr->find_f(drdmPair.first, 20);
			Lower.at(i) = Itr->find_f(drdmPair.second, 20);
		}
		return std::make_pair(std::move(Lower), std::move(Upper));
	}

	namespace detail {
		struct time_step {
		private:
			using this_type = time_step;
		private:
			unsigned int time;
			double f;
			double rm;
			double dr;
			double dm;
			bool Fail;
		public:
			time_step(unsigned int time_, double f_,double dr_, double dm_)
				:time(time_),f(f_),dr(dr_),dm(dm_), Fail(false){
				rm = dr / (dm + std::numeric_limits<double>::min());
			}
			void set(double f_, double dr_, double dm_) {
				f = f_;
				dr = dr_;
				dm = dm_;
				rm = dr / (dm + std::numeric_limits<double>::min());
			}
			void fail() { Fail = true; }
			void clear_fail() { Fail = false; }
			bool is_fail()const { return Fail; }
			unsigned int pos()const { return time; }
			double fraction()const { return f; }
			double reward()const { return f*dr; }
			double mortality()const { return f*dm; }
			double drdm()const { return dr / (dm + std::numeric_limits<double>::min()); }
			double lower_drdm()const { return (Fail ||f <= 0)? std::numeric_limits<double>::max(): dr / (dm + std::numeric_limits<double>::min()); }
			double upper_drdm()const { return (Fail ||f >= 1)? 0.: dr / (dm + std::numeric_limits<double>::min()); }
			friend bool operator<(const this_type& v1, const this_type& v2) {
				return v1.drdm() < v2.drdm();
			}
		};
	}
	template<typename info_iterator, typename fitness>
	state f_step_search_optimize(info_iterator InfoBeg, info_iterator InfoEnd, fitness&& Fitness, double Error) {
		unsigned int Size = std::distance(InfoBeg, InfoEnd);

		std::vector<double> Dif(Size, 0.1);

		std::vector<detail::time_step> TimeSet;
		auto Itr = InfoBeg;
		for (unsigned int i = 0; i<Size; ++i, ++Itr) {
			double f = 0.5;
			TimeSet.emplace_back(i, f, Itr->prey_reward(f), Itr->prey_mortality(f));
		}

		unsigned int ppos = std::numeric_limits<unsigned int>::max();
		bool Increase = true;
		while (std::any_of(TimeSet.begin(), TimeSet.end(), [](const detail::time_step& v) {return !v.is_fail(); })) {
			double R = std::accumulate(TimeSet.begin(), TimeSet.end(), 0., [](double v, const detail::time_step& t) {return v + t.reward(); });
			double M = std::accumulate(TimeSet.begin(), TimeSet.end(), 0., [](double v, const detail::time_step& t) {return v + t.mortality(); });
			auto lRM = std::min_element(TimeSet.begin(), TimeSet.end(), [](const detail::time_step& v1, const detail::time_step& v2) {return v1.lower_drdm() < v2.lower_drdm(); });
			auto uRM = std::max_element(TimeSet.begin(), TimeSet.end(), [](const detail::time_step& v1, const detail::time_step& v2) {return v1.upper_drdm() < v2.upper_drdm(); });

			if (is_prey_fitness_increasing(Fitness, R, M, uRM->upper_drdm())) {
				auto& t = *uRM;
				if (t.fraction() <= 0.0) {
					break;
				}
				unsigned int pos = t.pos();
				double or = t.reward();
				double om = t.mortality();
				auto Info = std::next(InfoBeg, pos);

				auto func = [=, &Fitness, &Info](double f){
					double dr = Info->prey_reward(f);
					double dm = Info->prey_mortality(f);
					if(is_prey_fitness_increasing(Fitness, R - or +f*dr, M - om + f*dm, dr/(dm+std::numeric_limits<double>::min())))return 1.0;
					else return -1.0;
				};

				//fail to increase minimum step
				if (func(t.fraction() + Error) < 0) {
					t.fail();
					continue;
				}

				if (func(1.0) > 0.0) {
					t.set(1.0, Info->prey_reward(1.0), Info->prey_mortality(1.0));
					t.fail();
				} else {
					auto Ans = boost::math::tools::bisect(func, t.fraction(), 1.0, boost::math::tools::eps_tolerance<double>(16));
					double f = (Ans.first + Ans.second)/2.0;
					t.set(f, Info->prey_reward(f), Info->prey_mortality(f));
					t.fail();
				}

				for (auto Itr = TimeSet.begin(); Itr != TimeSet.end(); ++Itr) {
					if (Itr == uRM)continue;
					Itr->clear_fail();
				}
			} else {
				auto& t = *lRM;
				if (t.fraction() <= 0.0) {
					break;
				}
				unsigned int pos = t.pos();
				double or = t.reward();
				double om = t.mortality();
				auto Info = std::next(InfoBeg, pos);

				auto func = [=, &Fitness, &Info](double f) {
					double dr = Info->prey_reward(f);
					double dm = Info->prey_mortality(f);
					if (is_prey_fitness_increasing(Fitness, R - or +f*dr, M - om + f*dm, dr / (dm + std::numeric_limits<double>::min())))return 1.0;
					else return -1.0;
				};

				//fail to increase minimum step
				if (func(t.fraction() - Error) > 0) {
					t.fail();
					continue;
				}

				if (func(0.0) < 0.0) {
					t.set(0.0, Info->prey_reward(0.0), Info->prey_mortality(0.0));
					t.fail();
				} else {
					auto Ans = boost::math::tools::bisect(func, 0.0, t.fraction(),boost::math::tools::eps_tolerance<double>(16));
					double f = (Ans.first + Ans.second) / 2.0;
					t.set(f, Info->prey_reward(f), Info->prey_mortality(f));
					t.fail();
				}

				for (auto Itr = TimeSet.begin(); Itr != TimeSet.end(); ++Itr) {
					if (Itr == lRM)continue;
					Itr->clear_fail();
				}
			}
		}

		state State;
		for (const auto& v : TimeSet)State.push_back(v.fraction());

		return State;
	}

	namespace detail {
		struct mutate {
			template<typename urbg>
			void operator()(state& x, urbg&& RandEng) {
				unsigned int pos = std::uniform_int_distribution<unsigned int>(0, x.size()-1)(RandEng);

				x[pos] = std::max(0.,std::min(1.,std::normal_distribution<double>(x[pos], 0.02)(RandEng)));
			}
		};
	}
	template<typename info_iterator, typename fitness>
	state hill_climbing_search(info_iterator InfoBeg, info_iterator InfoEnd, fitness&& Fitness, unsigned int StableStep, unsigned int MaxStep) {
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
		detail::mutate Mutate;
		hmLib::optimize::hill_climbing_search(State, Evaluate, Mutate, hmLib::optimize::breakers::limited_const_stable_breaker<state, double>(StableStep, MaxStep),hmLib::random::default_engine());
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
		state operator()(void) {
			return optimize_by_small_step();
		}
		std::pair<state, state> optimize_by_bisect(void) {
			return drdm_bisect_optimize(Container.begin(), Container.end(), PreyFitness);
		}
		state optimize_by_small_step(void) {
			return f_step_search_optimize(Container.begin(), Container.end(), PreyFitness, 1e-10);
		}
		state optimize_by_hill_climbing(unsigned int StableStep, unsigned int MaxStep) {
			return hill_climbing_search(Container.begin(), Container.end(), PreyFitness, StableStep, MaxStep);
		}

		const_iterator begin()const { return std::cbegin(Container); }
		const_iterator end()const { return std::cend(Container); }
		double get_prey_fitness(const state& x){
			double R = 0;
			double M = 0;
			for (unsigned int i = 0; i < x.size(); ++i) {
				R += x[i]*Container[i].prey_reward(x[i]);
				M += x[i] * Container[i].prey_mortality(x[i]);
			}

			return PreyFitness(R, M);
		}
	};

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
	template<typename hunt_prob>
	struct search_probforage_predation {
	private:
		hunt_prob Prob;	//function of successful hunting for given relative speed
		double s;	//average search time when the prey density is 1 (full predation) and swim speed is 1
		double h;	//average handling time for predation a prey
	public:
		search_probforage_predation(hunt_prob Prob_, double s_, double h_)
			: Prob(std::move(Prob_))
			, s(s_)
			, h(h_) {
		}
		double operator()(double f, double v, double u) const {
			if (f <= 0 || v <= 0)return 0.;
			return 1 / (s / v / f / Prob(v / (u + std::numeric_limits<double>::min())) + h);
		}
	};
	struct linear_prey_reward {
	private:
		double k;
	public:
		linear_prey_reward(double k_) :k(k_) {}
		double operator()(double u) const { return k*u; }
	};
	struct constant_predator_cost {
	private:
		double k;
	public:
		constant_predator_cost(double k_) :k(k_) {}
		double operator()(double v) const { return k; }
	};
	struct exp_diff_fitness {
	private:
		double alpha;	//speed of saturation
	public:
		exp_diff_fitness(double alpha_) :alpha(alpha_) {}
		double operator()(double r, double m) const {
			return 1 - std::exp(-alpha*r) - m;
		}
	};
}
#
#endif
