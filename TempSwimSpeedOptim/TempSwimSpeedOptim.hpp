#ifndef TEMPSWIMSPEEDOPTIM_INC
#define TEMPSWIMSPEEDOPTIM_INC
#
#include<tuple>
#include<vector>
#include<algorithm>
#include<numeric>
#include<boost/math/tools/roots.hpp>
#include<boost/math/tools/minima.hpp>
namespace tempss{
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
			return val(r + drdm*1e-10, m + 1e-10) > val(r - drdm*1e-10, m - 1e-10);
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
		template<typename prey_reward, typename predator_cost>
		time_info(const predation& Predation_, const prey_reward& PreyReward, const predator_cost& PredatorCost, double v_, double u_, double base_mu_, double d_)
			: Predation(Predation_)
			, v(v_)
			, u(u_)
			, base_mu(base_mu_)
			, d(d_) {
			c = PredatorCost(v);
			r = PreyReward(u);
			Max_drdm = std::max(0., PreyReward(0) / (prey_mortality(0) + std::numeric_limits<double>::min()));
			Min_drdm = std::max(0., PreyReward(1) / (prey_mortality(1) + std::numeric_limits<double>::min()));
			f_thr = find_predation_threshold(Predation_, v, u, c);
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
				min_drdm(),max_drdm(),
				boost::math::tools::eps_tolerance<double>(Bits)
			);
			return (val.first + val.second) / 2.0;
		}
		double max_drdm()const { return Max_drdm; }
		double min_drdm()const { return Min_drdm; }
		double f_threshold()const { return f_thr; }
	};

	template<typename info_iterator,typename fitness>
	auto drdm_bisect_optimize(info_iterator InfoBeg, info_iterator InfoEnd, fitness&& Fitness) {
		double MinRM = std::numeric_limits<double>::max();
		double MaxRM = std::numeric_limits<double>::lowest();
		for (auto Itr = InfoBeg; Itr != InfoEnd; ++Itr) {
			MaxRM = std::max(MaxRM, Itr->max_drdm());
			MinRM = std::min(MinRM, Itr->min_drdm());
		}

		auto Func = [=, &PreyFitness=Fitness](double drdm) {
			double r = 0;
			double m = 0;
			double max_drdm = 0;
			for (auto Itr = InfoBeg; Itr != InfoEnd; ++Itr) {
				double f = Itr->find_f(drdm,10);
				double dr = Itr->prey_reward(f);
				double dm = Itr->prey_mortality(f);

				r += dr*f;
				m += dm*f;

				if (f < 1) {
					max_drdm = std::max(max_drdm, dr / (dm + std::numeric_limits<double>::min()));
				}
			}
			return PreyFitness(r, m);
		};

		auto drdmPair = boost::math::tools::brent_find_minima(Func, MinRM, MaxRM, 10);

		unsigned int Size = std::distance(InfoBeg, InfoEnd);
		std::vector<double> Lower(Size, 0.);
		std::vector<double> Upper(Size, 0.);

		auto Itr = InfoBeg;
		for (unsigned int i = 0; i<Size; ++i, ++Itr) {
			Lower.at(i) = Itr->find_f(drdmPair.first,10);
			Upper.at(i) = Itr->find_f(drdmPair.second, 10);
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
			double dr;
			double dm;
		public:
			time_step(unsigned int time_, double f_,double dr_, double dm_):time(time_),f(f_),dr(dr_),dm(dm_){}
			void set(double f_, double dr_, double dm_) {
				f = f_;
				dr = dr_;
				dm = dm_;
			}
			unsigned int pos()const { return time; }
			double fraction()const { return f; }
			double reward()const { return f*dr; }
			double mortality()const { return f*dm; }
			double drdm()const { return dr / (dm + std::numeric_limits<double>::min()); }
			double lower_drdm()const { return f <= 0? std::numeric_limits<double>::max(): dr / (dm + std::numeric_limits<double>::min()); }
			double upper_drdm()const { return f >= 1? 0.: dr / (dm + std::numeric_limits<double>::min()); }
			friend bool operator<(const this_type& v1, const this_type& v2) {
				return v1.drdm() < v2.drdm();
			}
		};
	}
	template<typename info_iterator, typename fitness>
	auto f_step_search_optimize(info_iterator InfoBeg, info_iterator InfoEnd, fitness&& Fitness, std::vector<double> Lower, std::vector<double> Upper, double Error) {
		unsigned int Size = std::distance(InfoBeg, InfoEnd);
		if (Lower.size() != Size  || Size != Upper.size())throw std::exception();
		std::vector<detail::time_step> TimeSet;
	
		auto Itr = InfoBeg;
		for (unsigned int i = 0; i<Size; ++i, ++Itr) {
			double f = Lower[i];
			TimeSet.emplace_back(i, f, Itr->prey_reward(f), Itr->prey_mortality(f));
		}

		while (!TimeSet.empty()) {
			double R = std::accumulate(TimeSet.begin(), TimeSet.end(), 0., [](double v, const detail::time_step& t) {return v + t.reward(); });
			double M = std::accumulate(TimeSet.begin(), TimeSet.end(), 0., [](double v, const detail::time_step& t) {return v + t.mortality(); });
			auto lRM = std::min_element(TimeSet.begin(), TimeSet.end(), [](const detail::time_step& v1, const detail::time_step& v2) {return v1.lower_drdm() < v2.lower_drdm(); });
			auto uRM = std::max_element(TimeSet.begin(), TimeSet.end(), [](const detail::time_step& v1, const detail::time_step& v2) {return v1.upper_drdm() < v2.upper_drdm(); });

			if (is_prey_fitness_increasing(Fitness, R, M, uRM->upper_drdm())) {
				double lf = uRM->fraction();
				unsigned int pos = uRM->pos();
				double f = (Lower[pos] + Upper[pos]) / 2.0;
				Lower[pos] = lf;

				auto Info = std::next(InfoBeg, pos);
				double r = Info->prey_reward(f);
				double m = Info->prey_mortality(f);
				R += r - uRM->reward();
				M += m - uRM->mortality();

				if (Upper[pos] - Lower[pos] < Error) {
					TimeSet.erase(uRM);
				} else {
					uRM->set(f, r, m);
				}
			} else {
				double uf = lRM->fraction();
				unsigned int pos = lRM->pos();
				double f = (Lower[pos] + Upper[pos]) / 2.0;
				Upper[pos] = uf;

				auto Info = std::next(InfoBeg, pos);
				double r = Info->prey_reward(f);
				double m = Info->prey_mortality(f);
				R += r - uRM->reward();
				M += m - uRM->mortality();

				if (Upper[pos] - Lower[pos] < Error) {
					TimeSet.erase(lRM);
				} else {
					lRM->set(f, r, m);
				}
			}
		}

		return std::make_pair(Lower, Upper);
	}

	template<typename predation_, typename prey_fitness_>
	struct general_system {
		using predation = predation_;
		using prey_fitness = prey_fitness_;
		using state = std::vector<double>;
		using this_time_info = time_info<predation>;
	public:
		using container = std::vector<this_time_info >;
		using iterator = typename container::iterator;
		using const_iterator = typename container::const_iterator;
	private:
		predation Predation;
		prey_fitness PreyFitness;
		container Container;
	public:
		template<typename prey_reward, typename predator_cost, typename iterator>
		general_system(predation Predation_, prey_fitness PreyFitness_, prey_reward&& PreyReward, predator_cost&& PredatorCost, iterator VBeg, iterator VEnd, iterator UBeg, iterator UEnd, double base_mu_, double d_)
			: Predation(std::move(Predation_))
			, PreyFitness(std::move(PreyFitness_)) {
			if (std::distance(VBeg, VEnd) != std::distance(UBeg, UEnd))throw std::exception();

			for (; VBeg != VEnd; ++VBeg, ++UBeg) {
				Container.emplace_back(Predation, PreyReward, PredatorCost, *VBeg, *UBeg, base_mu_, d_);
			}
		}
		auto operator()(void) {
			auto Ans = drdm_bisect_optimize(Container.begin(), Container.end(), PreyFitness);
			auto Ans2 = f_step_search_optimize(Container.begin(), Container.end(), PreyFitness, std::move(Ans.first), std::move(Ans.second), 1e-10);

			return Ans2;
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
