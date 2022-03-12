#ifndef RANDOM_GENERATORS_H
#define RANDOM_GENERATORS_H

#include <limits>

namespace random {

	using std::numeric_limits;
	using real = std::false_t;
	using integer = std::true_t;

	namespace details {
		class Integer_Core {
		protected:
			struct pcg32_random_t {
				uint64_t state;
				uint64_t inc;
				pcg32_random_t(uint64_t init=get_current_nanoseconds()) : state(init), inc(init){}
			} rng;
			uint32_t pcg32_random_r()
			{
				uint64_t oldstate = rng.state;
		    // Advance internal state
				rng.state = oldstate * 6364136223846793005ULL + (rng.inc|1);
		    // Calculate output function (XSH RR), uses old state for max ILP
				uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
				uint32_t rot = oldstate >> 59u;
				return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
			}
		};
	};

	template<class T = real>
	class uniform : protected details::Integer_Core {

		double min,max;

	public:
		uniform_real(const double min=0.0,const double max=1.0) : min(min), max(max){}

		inline double operator()(){
			return min + (this->pcg32_random_r() / ( numeric_limits<uint32_t>::max() / (max-min) ) ) ;
		}

	};

	template<class T = integer, bool replace = false>
	class uniform : protected details::Integer_Core{
	// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
	// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

		uint32_t min_bound,max_bound;

	public:
		uniform(uint32_t max_bound=1) : min_bound(1),max_bound(max_bound){}
		uniform(uint32_t min_bound,uint32_t max_bound) : min_bound(min_bound),max_bound(max_bound){}

		uint32_t operator()(){
			return this->pcg32_random_r()%this->max_bound+this->min_bound;
		}

	};
	
	template<class T = integer, bool replace>
	class uniform : protected Integer_Core{

		vector<size_t> indices;

		void remove_index(uint32_t i){
			this->indices[(size_t)i]=this->indices.back();
			this->indices.pop_back();
		}

	public:
		uniform(uint32_t max_bound=0){
			this->indices.resize(max_bound);
			iota(this->indices.begin(),this->indices.end(),0);
		}

		uniform(int32_t min_bound,int32_t max_bound){
			this->indices.resize(std::abs(max_bound-min_bound+1));
			iota(this->indices.begin(),this->indices.end(),min_bound);
		}

		uint32_t operator()(){
			auto index = (size_t)this->pcg32_random_r()%this->indices.size();
			auto res = this->indices[index];
			this->remove_index(index);
			return res;
		}

	};
};

#endif