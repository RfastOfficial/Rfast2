#ifndef RANDOM_GENERATORS_H
#define RANDOM_GENERATORS_H

#include <limits>
#include <type_traits>
#include <chrono>
#include <vector>
#include <algorithm>
#include <zigg/header>

namespace Random
{

	using std::numeric_limits;
	using real = std::false_type;
	using integer = std::true_type;
	using std::iota;
	using std::vector;
	static zigg::Ziggurat ziggurat;

	namespace internal
	{

		static inline long long int get_current_nanoseconds()
		{
			return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		}

		class Integer_Core
		{
		public:
			// The next 3 lines are the requirements for UniformRandomBitGenerator.
			using result_type = uint32_t;
			static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
			static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

		protected:
			struct pcg32_random_t
			{
				uint64_t state;
				uint64_t inc;
				pcg32_random_t(uint64_t init = get_current_nanoseconds()) : state(init), inc(init) {}
			} rng;
			result_type pcg32_random_r()
			{
				uint64_t oldstate = rng.state;
				// Advance internal state
				rng.state = oldstate * 6364136223846793005ULL + (rng.inc | 1);
				// Calculate output function (XSH RR), uses old state for max ILP
				result_type xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
				result_type rot = oldstate >> 59u;
				return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
			}
		};
	}

	template <class T, bool replace = false>
	class uniform : public internal::Integer_Core
	{

		vector<size_t> indices;

		void remove_index(result_type i)
		{
			this->indices[(size_t)i] = this->indices.back();
			this->indices.pop_back();
		}

	public:
		uniform(result_type max_bound = 0)
		{
			this->indices.resize(max_bound);
			iota(this->indices.begin(), this->indices.end(), 0);
		}

		uniform(int32_t min_bound, int32_t max_bound)
		{
			this->indices.resize(std::abs(max_bound - min_bound + 1));
			iota(this->indices.begin(), this->indices.end(), min_bound);
		}

		result_type operator()()
		{
			auto index = (size_t)this->pcg32_random_r() % this->indices.size();
			auto res = this->indices[index];
			this->remove_index(index);
			return res;
		}
	};

	template <class T>
	class uniform<T, true> : public internal::Integer_Core
	{
		// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
		// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

		result_type min_bound, max_bound;

	public:
		uniform(result_type max_bound = 1) : min_bound(1), max_bound(max_bound) {}
		uniform(result_type min_bound, result_type max_bound) : min_bound(min_bound), max_bound(max_bound) {}

		result_type operator()()
		{
			return this->pcg32_random_r() % this->max_bound + this->min_bound;
		}
	};

	template <>
	class uniform<real, false> : public internal::Integer_Core
	{

		const double min, max;

	public:
		uniform(const double min = 0.0, const double max = 1.0) : min(min), max(max) {}

		inline double operator()()
		{
			return min + (this->pcg32_random_r() * (max - min) / internal::Integer_Core::max());
		}
	};

	static uniform<real> rng(0, 1);

	class Gamma
	{
		double rate, d, c, boosted_shape, inv_init_shape;
		bool boosted;

	public:
		Gamma(double shape, double rate) : rate(1.0 / rate)
		{
			if (shape < 1.0)
			{
				boosted_shape = shape + 1.0;
				boosted = true;
				inv_init_shape = 1.0 / shape; // Precompute to avoid division
			}
			else
			{
				boosted_shape = shape;
				boosted = false;
			}
			d = boosted_shape - 1.0 / 3.0;
			c = 1.0 / std::sqrt(9.0 * d);
		}

		double operator()()
		{
			while (true)
			{
				double x = ziggurat.rnorm();
				double v = 1.0 + c * x;
				v = v * v * v;
				if (v > 0)
				{
					double u = rng();
					double x2 = x * x;
					if (u < 1.0 - 0.0331 * x2 * x2 || std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v)))
					{
						double res = d * v * rate;
						if (boosted)
						{
							return res * std::exp(std::log(u) * inv_init_shape); // Faster than pow()
						}
						else
						{
							return res;
						}
					}
				}
			}
		}
	};

	// class Gamma
	// {
	// 	double shape, rate, d, c;

	// public:
	// 	Gamma(double shape, double rate) : shape(shape), rate(1.0 / rate), d(shape - 1.0 / 3.0), c(1.0 / std::sqrt(9.0 * d)) {}

	// 	double operator()()
	// 	{
	// 		// Marsaglia and Tsang method for rate >= 1
	// 		while (true)
	// 		{
	// 			double x = ziggurat.rnorm(), x2 = x * x;
	// 			double v = 1.0 + c * x;
	// 			v = v * v * v;
	// 			double u = rng();

	// 			if (v > 0 && (u < 1.0 - 0.0331 * x2 * x2 || std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v))))
	// 			{
	// 				return d * v * rate;
	// 			}
	// 		}
	// 	}
	// };

	class BetaOne : public Gamma
	{

	public:
		BetaOne(double alpha) : Gamma(alpha, 1.0) {}

		inline double operator()()
		{
			double x = Gamma::operator()();
			return x / (x + Gamma::operator()());
		}
	};

	class Beta : public BetaOne
	{
		Gamma beta_d;

	public:
		Beta(double alpha, double beta) : BetaOne(alpha), beta_d(Gamma(beta, 1.0)) {}

		inline double operator()()
		{
			double x = Gamma::operator()();
			return x / (x + beta_d());
		}
	};

	class Exp : public Gamma
	{
	public:
		Exp(double rate) : Gamma(1.0, rate) {}
	};

	class Chisq : public Gamma
	{
	public:
		// rate must be 2 but Gamma divides with 1. So to undo it we need to divided first with 1 and pass it.
		Chisq(double df) : Gamma(df / 2.0, 1.0 / 2.0) {}
	};

	class Geom : public uniform<real, false>
	{
		double lambda;

	public:
		Geom(double prob) : lambda(-log(1 - prob)) {}

		inline double operator()()
		{
			return std::floor(std::log(uniform<real, false>::operator()()) / lambda);
		}
	};

	class Cauchy
	{
		double location, scale;

	public:
		Cauchy(double location = 0, double scale = 1) : location(location), scale(scale) {}

		double operator()()
		{
			return location + scale * std::tan(M_PI * ziggurat.rnorm());
		}
	};

	class StudentT : public Chisq
	{
		double sqrt_df, ncp_sqrt_df;

	public:
		StudentT(double df, double ncp = 0) : Chisq(df), sqrt_df(std::sqrt(df)), ncp_sqrt_df(ncp * sqrt_df) {}

		double operator()()
		{
			// return (ziggurat.rnorm() + ncp) / std::sqrt(Chisq::operator()() / df);
			return (ziggurat.rnorm() * sqrt_df + ncp_sqrt_df) / std::sqrt(Chisq::operator()());
		}
	};
}

#endif