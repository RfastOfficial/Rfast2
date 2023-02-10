/// [[Rcpp::depends(RcppArmadillo)]];
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <chrono>
#include <string>
#include "templates.h"
#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace Rcpp;

#include <windows.h>
#include <psapi.h>

class PC {
    ULARGE_INTEGER lastCPU, lastSysCPU, lastUserCPU;
    int numProcessors;
    HANDLE self;

public:
    PC(){
        SYSTEM_INFO sysInfo;
        FILETIME ftime, fsys, fuser;

        GetSystemInfo(&sysInfo);
        numProcessors = sysInfo.dwNumberOfProcessors;

        GetSystemTimeAsFileTime(&ftime);
        memcpy(&lastCPU, &ftime, sizeof(FILETIME));

        self = GetCurrentProcess();
        GetProcessTimes(self, &ftime, &ftime, &fsys, &fuser);
        memcpy(&lastSysCPU, &fsys, sizeof(FILETIME));
        memcpy(&lastUserCPU, &fuser, sizeof(FILETIME));
    }
    double getCurrentCPU(){
        FILETIME ftime, fsys, fuser;
        ULARGE_INTEGER now, sys, user;
        double percent;

        GetSystemTimeAsFileTime(&ftime);
        memcpy(&now, &ftime, sizeof(FILETIME));

        GetProcessTimes(self, &ftime, &ftime, &fsys, &fuser);
        memcpy(&sys, &fsys, sizeof(FILETIME));
        memcpy(&user, &fuser, sizeof(FILETIME));
        percent = (sys.QuadPart - lastSysCPU.QuadPart) +
        (user.QuadPart - lastUserCPU.QuadPart);
        percent /= (now.QuadPart - lastCPU.QuadPart);
        percent /= numProcessors;
        lastCPU = now;
        lastUserCPU = user;
        lastSysCPU = sys;

        return percent * 100;
    }
    SIZE_T getCurrentVirtualMemory(){
        return getMemoryInfo().PrivateUsage;
    }
    SIZE_T getCurrentPhysicalMemory(){
        return getMemoryInfo().WorkingSetSize;
    }

private:

    PROCESS_MEMORY_COUNTERS_EX getMemoryInfo(){
        PROCESS_MEMORY_COUNTERS_EX pmc;
        GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
        return pmc;
    }
};

SIZE_T distance(SIZE_T a,SIZE_T b){
    return a > b ? a-b : b-a;
}

static NumericVector measure_time(SEXP expr,SEXP env,const int tim,PC &pc){
    steady_clock sc;
    NumericVector times(tim);
    SIZE_T memory = 0;
    double cpu = 0.0;
    double sum_t=0,max_t,min_t;
    for (int i = 0; i < tim; ++i){
        SIZE_T mem = pc.getCurrentPhysicalMemory();
        auto start = sc.now();
        Rf_eval(expr,env);
        auto end = sc.now();
        memory += distance(pc.getCurrentPhysicalMemory(),mem);
        cpu += pc.getCurrentCPU();
        times[i] = static_cast<chrono::duration<double>>(end - start).count();
        sum_t+=times[i];
    }
    min_max<double>(&times[0],&times[tim-1]+1,min_t,max_t);
    return NumericVector::create(min_t,sum_t/tim,max_t,cpu/tim,memory/tim);
}

//[[Rcpp::export]]
NumericMatrix benchmark(List exprs,SEXP env,const int tim,IntegerVector indices){
    NumericMatrix res(exprs.length(),5);
    PC pc;
    for(auto& index : indices){
        res.row(index-1) = measure_time(exprs[index-1],env,tim,pc);
    }
    return res;
}

RcppExport SEXP Rfast2_benchmark(SEXP exprsSEXP,SEXP envSEXP,SEXP timSEXP,SEXP indicesSEXP){
    BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< List  >::type exprs(exprsSEXP);
    traits::input_parameter< SEXP  >::type env(envSEXP);
    traits::input_parameter< const int  >::type tim(timSEXP);
    traits::input_parameter< IntegerVector  >::type indices(indicesSEXP);
    __result = benchmark(exprs,env,tim,indices);
    return __result;
    END_RCPP
}








/*
using namespace chrono;

//[[Rcpp::export]]
string tic_c(){
constexpr auto zero_d = time_point<steady_clock>(seconds(0));
double s = static_cast<duration<double>>(steady_clock().now()-zero_d).count();
return to_string(s);
}

//[[Rcpp::export]]
NumericVector duration_s(string x,string y){
double d=stod(x)-stod(y);
CharacterVector unit="seconds";
NumericVector f(1);
if(d >= 1.0){ // seconds
unit = "seconds";
}
else if((d*=1e+03) >= 1.0){ //milliseconds
unit = "milliseconds";
}else{//microseconds
unit = "microseconds";
d *= 1e+03;
}
f.names()=unit;
f[0]=d;
return f;
}

*/
