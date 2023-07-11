#include "Rfast2/templates.h"

using namespace std;

#ifndef FUNCS_TEMPLATES_H
#define FUNCS_TEMPLATES_H

/*
* sum(Bin_Fun( a, b))
* ? is a binary function
* a, b are iterators (should be same size)
*/
template<Binary_Function B, class T1, class T2>
inline double apply_funcs(T1 a, T2 b, const int sz){
	double ret = 0.0;
	for(int i=0;i<sz;i++,++a,++b){
		ret+=B(*a, *b);
	}
	return ret;
}

/*
* sum(fun0( Bin_Fun( fun1(a), fun2(b))))
* fun0, fun1, fun2 are Unary_Functions
* ? is a binary function
* a, b are iterators (should be same size)
*/
template<Unary_Function F0, Unary_Function F1, Binary_Function B, Unary_Function F2, typename T1, typename T2>
inline double apply_funcs(T1 a, T2 b, const int sz){
	double ret = 0.0;
	for(int i=0;i<sz;i++,++a,b++){
		ret+=F0( B(F1(*a), F2(*b)) );
	}
	return ret;
}

/*
* sum(fun0( fun1(a) ? fun2(b))) (only now, one of fun0, fun1, fun2 doesn't exist)
* fun0, fun1, fun2 are Unary_Functions
* ? is a binary function
* a, b are iterators (should be same size)
* type indicates which of fun0, fun1, fun2 doesn't exist
*/
template<Unary_Function F0, Unary_Function F1, Binary_Function B, typename T1, typename T2>
inline double apply_funcs(T1 a, T2 b, const int sz, const int type){
	double ret = 0.0;
	if(type==0){
		for(int i=0;i<sz;i++,++a,++b){
			ret+=B(F0(*a), F1(*b));
		}
	}
	else if(type==1){
		for(int i=0;i<sz;i++,++a,++b){
			ret+=F0( B(*a, F1(*b)) );
		}
	}
	else{
		for(int i=0;i<sz;i++,++a,++b){
			ret+=F0( B(F1(*a), *b) );
		}
	}
	return ret;
}

/*
* sum(fun0( fun1(a) ? fun2(b))) (like before, only now only one fun exists)
* fun0, fun1, fun2 are Unary_Functions
* ? is a binary function
* a, b are iterators (should be same size)
* type indicates which of fun0, fun1, fun2 exists
*/
template<Unary_Function F, Binary_Function B, typename T1, typename T2>
inline double apply_funcs(T1 a, T2 b, const int sz, const int type){
	double ret = 0.0;
	if(type==0){
		for(int i=0;i<sz;i++,++a,++b){
			ret+=F( B(*a, *b) );
		}
	}
	else if(type==1){
		for(int i=0;i<sz;i++,++a,++b){
			ret+=B(F(*a), *b);
		}
	}
	else{
		for(int i=0;i<sz;i++,++a,++b){
			ret+=B(*a, F(*b));
		}
	}
	return ret;
}

template<Binary_Function B, typename T1, typename T2, typename T3>
inline void apply_funcs(T1 a, T2 b, T3 out, const int sz){

  for(int i=0;i<sz;i++,++a,++b,++out){
    *out=B(*a, *b);
  }
  return;
}

/*
* each element of out is filled with (fun0( fun1(a) ? fun2(b)))
* fun0, fun1, fun2 are Unary_Functions
* ? is a binary function
* a, b, out are iterators (should refer to elements with the same size)
*/
template<Unary_Function F0, Unary_Function F1, Binary_Function B, Unary_Function F2, typename T1, typename T2, typename T3>
inline void apply_funcs(T1 a, T2 b, T3 out, const int sz){
	for(int i=0;i<sz;i++,++a,++b, ++out){
		*out=F0( B(F1(*a), F2(*b)) );
	}
	return;
}

/*
*  each element of out is filled with sum(fun0( fun1(a) ? fun2(b))) (only now, one of fun0, fun1, fun2 doesn't exist)
* fun0, fun1, fun2 are Unary_Functions
* ? is a binary function
* a, b are iterators (should be same size)
* type indicates which of fun0, fun1, fun2 doesn't exist
*/
template<Unary_Function F0, Unary_Function F1, Binary_Function B, typename T1, typename T2, typename T3>
inline void apply_funcs(T1 a, T2 b, T3 out, const int sz, const int type){
	if(type==0){
		for(int i=0;i<sz;i++,++a,++b,++out){
			*out=B(F0(*a), F1(*b));
		}
	}
	else if(type==1){
		for(int i=0;i<sz;i++,++a,++b,++out){
			*out=F0( B(*a, F1(*b)) );
		}
	}
	else{
		for(int i=0;i<sz;i++,++a,++b,++out){
			*out=F0( B(F1(*a), *b) );
		}
	}
	return;
}

/*
*  each element of out is filled with sum(fun0( fun1(a) ? fun2(b))) (like before, only now only one fun exists)
* fun0, fun1, fun2 are Unary_Functions
* ? is a binary function
* a, b are iterators (should be same size)
* type indicates which of fun0, fun1, fun2 exists
*/
template<Unary_Function F, Binary_Function B, typename T1, typename T2, typename T3>
inline void apply_funcs(T1 a, T2 b, T3 out, const int sz, const int type){
	if(type==0){
		for(int i=0;i<sz;i++,++a,++b,++out){
			*out=F( B(*a, *b) );
		}
	}
	else if(type==1){
		for(int i=0;i<sz;i++,++a,++b,++out){
			*out=B(F(*a), *b);
		}
	}
	else{
		for(int i=0;i<sz;i++,++a,++b,++out){
			*out=B(*a, F(*b));
		}
	}
	return;
}

#endif
