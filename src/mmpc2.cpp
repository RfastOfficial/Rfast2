
// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include "reg_lib2.h"
#include <chrono>
#include "system_files.h"

using namespace std;
using namespace arma;

#define cbind join_rows
#define rbind join_cols

template<class Ret,class T>
void combn(T vals, const int n, const unsigned int start_idx,
           std::vector<double>& combn_data, Ret& combn_ds, unsigned int& combn_col) {
    if (!n) {
        for (unsigned int i = 0; i < combn_ds.n_rows && combn_col < combn_ds.n_cols; i++) {
            combn_ds(i, combn_col) = combn_data[i];
        }
        ++combn_col;
        return;
    }
    for (unsigned int i = start_idx; i <= (vals.size() - n); i++) {
        combn_data[combn_ds.n_rows - n] = vals[i];
        combn<Ret,T>(vals, n - 1, i + 1, combn_data, combn_ds, combn_col);
    }
}

template<class Ret,class T>
Ret find_combn(T vals, const int n) {
    static unsigned int combn_col = 0;
    const unsigned int nrows = n;
    const unsigned int ncols = std::round(R::choose(vals.size(), n));
    Ret combn_ds(nrows, ncols);
    std::vector<double> combn_data(nrows);
    const unsigned int start_idx = 0;
    combn_col = 0;
    combn<Ret,T>(vals, n, start_idx, combn_data, combn_ds, combn_col);
    return combn_ds;
}

template<class T>
T pmax(T x,T y){
    for(unsigned int i=0;i<x.n_elem;++i)
        x[i]=max(x[i],y[i]);
    return x;
}

//[[Rcpp::export]]
arma::mat quasi_poisson_only(arma::mat x, arma::colvec y, const double ylogy, const double tol,const int maxiters){
    const unsigned int n=x.n_rows,pcols=x.n_cols,d=2;
    unsigned int i;
    int ij;

    colvec b_old(d),b_new(d),L1(d),yhat(n);
    mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones);
    vec m(n),z_col_1(n);
    mat F(2,pcols);
    double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(mean(y));
    for(i=0;i<pcols;++i){
        b_old(0)=lgmeany;
        b_old(1)=0;
        z_col_1=x.col(i);
        z.col(1)=z_col_1;
        z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
        ij=2;
        for(dif=1.0;dif>tol;){
            sm=szm=sz2m=0.0;
            yhat=z*b_old;
            m=(exp(yhat));
            L1=z_tr*(y-m);
            sm=sum(m);
            szm=sum(m%z_col_1);
            sz2m=sum(m%arma::square(z_col_1));
            t=1.0/(sm*sz2m-szm*szm);
            inv_L2.at(0,0)=sz2m*t;
            inv_L2.at(0,1)=inv_L2.at(1,0)=-szm*t;
            inv_L2.at(1,1)=sm*t;
            b_new=b_old+inv_L2*L1;
            dif=sum(abs(b_new-b_old));
            b_old=b_new;
            if(++ij==maxiters)
                break;
        }
        F(0,i)= 2.0*(ylogy-sum(y%yhat));
        F(1,i) = sum(arma::square(y-m)/m)/(n-d);
    }
    return F;
}
colvec poisson_only(mat x, colvec y,const double ylogy){
    const unsigned int n=x.n_rows,pcols=x.n_cols,d=2;
    unsigned int i;
    colvec b_old(d),b_new(d),L1(d),yhat(n);
    mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones);
    vec m(n),z_col_1(n);
    NumericVector F(pcols);
    double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(mean(y));
    for(i=0;i<pcols;++i){
        b_old[0]=lgmeany;
        b_old[1]=0;
        z_col_1=x.col(i);
        z.col(1)=z_col_1;
        z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
        for(dif=1.0;dif>0.000000001;){
            sm=szm=sz2m=0.0;
            yhat=z*b_old;
            m=exp(yhat);
            L1=z_tr*(y-m);
            sm=sum(m);
            szm=sum(m%z_col_1);
            sz2m=sum(m%square(z_col_1));
            t=1.0/(sm*sz2m-szm*szm);
            inv_L2.at(0,0)=sz2m*t;
            inv_L2.at(0,1)=inv_L2.at(1,0)=-szm*t;
            inv_L2.at(1,1)=sm*t;
            b_new=b_old+inv_L2*L1;
            dif=sum(abs(b_new-b_old));
            b_old=b_new;
        }
        F[i]=2.0*(ylogy-sum(y%yhat));
    }
    return F;
}

colvec logistic_only(mat x, colvec y,const double d0,const double my,const double sy,const double lgmy,const double tol){
    int maxiters = 100;
    const unsigned int n=x.n_rows,pcols=x.n_cols;

    unsigned int i=0;
    int j=0;
    colvec be(2),expyhat(n),W(n,fill::zeros),x_col(n),x2_col(n),de(n);
    mat yhat(n,1),p(n,1);
    colvec F(pcols);
    double d1=0,d2=0,t=0,dera2=0.0,sp=0.0,derb=0.0,dera=0.0,derab=0.0,derb2=0.0;

    double W0=my*(1-my);
    double dera20=n*W0;
    colvec de0(n);
    de0=y-my;
    double dera0=0;

    for(i=0;i<pcols;++i){
        d1=d0;
        be[0]=lgmy;
        be[1]=0;
        x_col=x.col(i);
        x2_col=square(x_col);
        derb=sum(de0%x_col);
        derab=W0*sum(x_col);
        derb2=W0*sum(x2_col);
        t=dera20 * derb2 - derab*derab;
        be[0]=be[0]+(derb2 * dera0 - derab * derb)/t;
        be[1]=be[1]+( - derab * dera0 + dera20 * derb )/t;

        yhat.col(0) = be[0]+be[1]*x_col;
        expyhat=exp(-yhat.col(0));
        p.col(0) = 1 / ( 1 + expyhat );

        d2 = -2*calcDevRes(p,y,yhat);
        //errors("d1",d1);
        //errors("d2",d2);
        j=2;
        while(j++<maxiters && d1-d2>tol){
            d1=d2;
            W=p%(1-p);
            dera2=sum(W);
            sp=sum(p.col(0));
            de=y-p.col(0);
            dera=sy-sp;
            derb=sum(de%x_col);
            derab=sum(W%x_col);
            derb2=sum(W%x2_col);
            t=dera2 * derb2 - derab*derab;
            be[0]=be[0]+(derb2 * dera - derab * derb)/t;
            be[1]=be[1]+( - derab * dera + dera2 * derb )/t;
            yhat.col(0) = be[0]+be[1]*x_col;
            expyhat=exp(-yhat.col(0));
            p.col(0) = 1 / ( 1 + expyhat );

            d2 = -2*calcDevRes(p,y,yhat);
        }

        F(i) = d2;
    }
    return F;
}

uvec std_setdiff(uvec& x,uvec b) {
    //panta x,y unique
    vector<unsigned int> out;
    std::set_difference(x.begin(), x.end(), b.begin(), b.end(),back_inserter(out));

    return conv_to<uvec>::from(out);
}

uvec std_setdiff(uvec& x,unsigned int y) {
    //panta x,y unique
    vector<unsigned int> out,b{y};
    std::set_difference(x.begin(), x.end(), b.begin(), b.end(),back_inserter(out));

    return conv_to<uvec>::from(out);
}

ucolvec find_arr_indices_c(umat x){
    uvec inds = find(x);
    const double nr=(double)x.n_rows;
    ucolvec res(inds.n_elem);
    for(unsigned int i=0;i<inds.n_elem;++i){
        res(i) = floor(inds(i)/nr);
    }
    return res;
}


double calc_devi_0(vec y,mat x,const add_term_ini_vars& ini,const int maxiters,const double tol,string test,const double lgmy=0.0,const double my=0.0,const double sy=0.0,const double d1=0.0,const double ylogy=0.0){
    double devi_0;
    if (test == "poisson" || test == "qpoisson") {  //## Poisson regression
        devi_0=2*ylogy+glm_poisson3(x, y,lgmy,tol,maxiters);
    } else if(test == "logistic") {
        devi_0=glm_logistic3(x,y,&ini.my[0],ini.ini,tol,maxiters);
    }else{
        stop("Error: wrong type.\n");
    }
    return devi_0;
}

mat Cbind(colvec& oness,mat x){
    x.insert_cols(0,oness);
    return x;
}

void push_back(uvec &x,unsigned int val){
    uvec y(1);
    y(0)=val;
    x.insert_rows(x.n_elem,y.row(0));
}


mat rep(unsigned int nrow,unsigned int ncol,unsigned int val){
    mat f(nrow,ncol);
    f.fill(val);
    return f;
}

List mmpc2(vec y,mat x,int max_k = 3,const double threshold = 0.05,const string test = "logistic",SEXP Ini=R_NilValue,const bool parallel = false,const int maxiters = 100,const double tol = 1e-07,const bool backward = false) {

    Timer timer;

    timer.Start();
    add_term_ini_vars ini = add_term_ini(y, test, tol, maxiters);
    int sel,dm;
    umat cand;
    double devi_0,m;
    double my=0.0,sy=0.0,lgmy=0.0,d1=0.0,ylogy=0.0;

    colvec pval,oness=ones(x.n_rows);
    const int varsize = x.n_cols;

    if ( max_k < 1L )
        stop("invalid max_k option");
    if ( max_k > varsize )
        max_k = varsize;
    if ( threshold <= 0 || threshold > 1 )
        stop("invalid threshold option");

    unsigned int ntests = 0;
    List kapa_pval(max_k);
    double alpha = log(threshold);

    if ( Rf_isNull(Ini) ) {
        //errors("Enter null ini, 259\n");
        colvec stat,mod;
        const int n=x.n_rows;
        pval=colvec(varsize);
        if (test == "logistic") {
            //errors("Logistic test 264\n");
            my=mean(y);
            sy=my*n;
            lgmy=log(my/(1-my));
            d1 = -2*(sy*log(my)+(n-sy)*log(1-my));
            //errors(d1);
            mod = logistic_only(x, y,d1,my,sy,lgmy,tol);
            //mod.print("mod");
            stat = d1 - mod;
            for(int i=0;i<varsize;++i)
                pval[i] = R::pchisq(stat[i], 1, 0,1);
            //errors("logistic end\n");
        }
        else if (test == "poisson") {
            m = mean(y);
            lgmy=log(m);
            ylogy=calcylogy(y,n);
            d1 = 2 * ylogy - 2 * n * m * lgmy;
            mod = poisson_only(x, y,ylogy);
            stat = d1 - mod;
            for(int i=0;i<varsize;++i)
                pval[i] = R::pchisq(stat[i], 1,0,1);
        }
        else if (test == "qpoisson") {
            m = mean(y);
            lgmy=log(m);
            ylogy=calcylogy(y,n);
            d1 = 2 * ylogy - 2 * n * m * lgmy;
            mat mod = quasi_poisson_only(x, y, ylogy, tol,maxiters);
            rowvec stat = (d1 - mod.row(0))/mod.row(1);
            for(int i=0;i<varsize;++i)
                pval[i] = R::pf(stat[i], 1, n - 2,0,1);
        }

        ntests = varsize;
    } else {
        pval = as<colvec>(Ini);
    }
    //errors(300);
    colvec univ=pval;
    uvec vars = find(pval < alpha);
    uvec sela(1);
    if ( vars.n_elem > 0 ) {
        sela(0)=pval.index_min();
    } else
        sela=vars;
    //errors(308);
    uvec ide;
    vars = std_setdiff(vars, sela);
    //errors(311);
    colvec pval2;
    ////## 1 selected
    if ( vars.n_elem > 0  &&  max_k >= 1 ) {
        //errors("Enter 1 selected");
        StringVector a(max_k);
        for(int i=0;i<max_k;++i)
            a[i]="kappa="+to_string(i+1);
        //kapa_pval.names()=a;
        devi_0 = calc_devi_0(y,Cbind(oness,x.col(sela(0))),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
        //errors(336," ",devi_0);
        pval2 = add_term_c(y, Cbind(oness,x.col(sela(0))),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
        //errors(339);
        //errors("vars.n_elem: ",vars.n_elem);
        /*sela.print("sela");
        vars.print("vars");
        pval2.print("pval2");*/
        kapa_pval[0] = join_rows((mat)pval2, join_rows(conv_to<colvec>::from(vars), (mat)repmat(conv_to<mat>::from(sela),vars.n_elem,1))+1);
        //errors(340);
        ntests += vars.n_elem;
        pval(vars) = pmax<colvec>(pval(vars), pval2);
        ide = find(pval(vars) < alpha);
        vars = vars(ide);
        if(vars.n_elem){
            sel = pval(vars).index_min();
            //errors("push_back");
            push_back(sela,vars(sel));
            //errors("aaaaaaaaaaaaa ",vars(sel));
            vars = std_setdiff(vars, vars(sel));
        }
        //errors("end 1 selected\n");
    }  //////## end  if ( vars.n_elem > 0  &  max_k >= 1 ) {

    ////## 2 selected
    if ( vars.n_elem > 0  &&  max_k >= 1 ) {
        //errors("Enter 2 selected\n");
        devi_0 = calc_devi_0(y,Cbind(oness,x.col(sela(1))),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
        pval2 = add_term_c(y, Cbind(oness,x.col(sela(1))),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
        //errors(369);
        kapa_pval[0] = join_cols( as<mat>(kapa_pval[0]), join_rows((mat)pval2, join_rows(conv_to<colvec>::from(vars), rep(vars.n_elem,1,sela(1)))+1));
        //errors(371);
        ntests +=vars.n_elem;
        pval(vars) = pmax<colvec>(pval(vars), pval2);
        ide = find(pval(vars) < alpha);
        vars = vars(ide);
        if ( vars.n_elem > 0  &&  max_k >= 2 ) {
            devi_0 = calc_devi_0(y,Cbind(oness,x.cols(sela)),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
            pval2 = add_term_c(y, Cbind(oness,x.cols(sela)),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
            //errors(389);
            kapa_pval[1] =join_rows( (mat)pval2, join_rows(conv_to<colvec>::from(vars), (mat)repmat(conv_to<rowvec>::from(sela), vars.n_elem,1))+1);
            //errors(390);
            pval(vars)=pmax<colvec>(pval(vars),pval2);
            ide = find(pval(vars) < alpha);
            ide=find(pval(vars)<alpha);
            vars = vars(ide);
        }  //## end  if ( max_k >= 2 ) {
        if( vars.n_elem){
            sel=pval(vars).index_min();
            push_back(sela,vars(sel));
            vars = std_setdiff(vars, vars(sel));
            dm = sela.n_elem;
        }
        //errors("end 2 selected\n");
    }  ////## end  if ( vars.n_elem > 0  &  max_k >= 1 ) {

    //## 3 selected
    if ( vars.n_elem > 0  &&  max_k >= 1 ) {
        //errors("Enter 3 selected\n");
        if ( max_k >= 1 ) {
            devi_0 = calc_devi_0(y,Cbind(oness,x.col(sela(2))),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
            pval2 = add_term_c(y, Cbind(oness,x.col(sela(2))),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
            kapa_pval[0] = rbind( as<mat>(kapa_pval[0]), cbind(pval2,cbind(conv_to<mat>::from(vars), rep(vars.n_elem,1,sela(2)))+1));
            //errors(411);
            ntests+=vars.n_elem;
            pval(vars)=pmax<colvec>(pval(vars),pval2);
            ide = find(pval(vars) < alpha);
            ide=find(pval(vars)<alpha);
            vars = vars(ide);
        }  //## end  if ( max_k >= 1 ) {
        if ( vars.n_elem > 0  &&  max_k >= 2 ) {
            cand = find_combn<umat,uvec>(sela, 2); // w
            devi_0 = calc_devi_0(y,Cbind(oness,x.cols(cand.col(1))),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
            pval2 = add_term_c(y, Cbind(oness,x.cols(cand.col(1))),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
            //errors(423);
            kapa_pval[1] = join_cols( as<mat>(kapa_pval[1]), join_rows(pval2,join_rows(conv_to<mat>::from(vars), (mat)repmat(conv_to<mat>::from(cand.col(1).t()), vars.n_elem,1))+1));
            //errors(427);
            ntests+=vars.n_elem;
            pval(vars)=pmax<colvec>(pval(vars),pval2);
            ide = find(pval(vars) < alpha);
            ide=find(pval(vars)<alpha);
            vars = vars(ide);
            //errors(431);
            devi_0 = calc_devi_0(y,Cbind(oness,x.cols(cand.col(2))),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
            pval2 = add_term_c(y, Cbind(oness,x.cols(cand.col(2))),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
            //errors(433);
            kapa_pval[1] = join_cols( as<mat>(kapa_pval[1]), join_rows( pval2, join_rows(conv_to<mat>::from(vars), (mat)repmat(conv_to<mat>::from(cand.col(2).t()), vars.n_elem,1))+1));
            //errors(436);
            ntests+=vars.n_elem;
            pval(vars)=pmax<colvec>(pval(vars),pval2);
            ide = find(pval(vars) < alpha);
            ide=find(pval(vars)<alpha);
            vars = vars(ide);
        }  //## end  if ( max_k >= 2 ) {
        if ( vars.n_elem > 0  &&  max_k >= 3 ) {
            //errors(444);
            devi_0 = calc_devi_0(y,Cbind(oness,x.cols(sela)),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
            pval2 = add_term_c(y, Cbind(oness,x.cols(sela)),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
            //errors(447);
            kapa_pval[2] = join_rows( (mat)pval2, join_rows(conv_to<mat>::from(vars), (mat)repmat(conv_to<rowvec>::from(sela), vars.n_elem,1))+1);
            ntests+=vars.n_elem;
            pval(vars)=pmax<colvec>(pval(vars),pval2);
            ide = find(pval(vars) < alpha);
            ide=find(pval(vars)<alpha);
            vars = vars(ide);
        }
        //## end  if ( max_k >= 3 ) {
        if( vars.n_elem > 0){
            sel=pval(vars).index_min();
            push_back(sela,vars(sel));
            vars = std_setdiff(vars, vars(sel));
        }
        //errors("end 3 selected\n");
    }  //## end  if ( vars.n_elem > 0  &  max_k >= 2 ) {

    //## 4 selected
    while ( vars.n_elem > 0  &&  max_k >= 1 ) {
        dm = sela.n_elem-1;
        devi_0 = calc_devi_0(y,Cbind(oness,x.col(sela(dm))),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
        pval2 = add_term_c(y, Cbind(oness,x.col(sela(dm))),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
        //errors(468);
        kapa_pval[0] = join_cols( as<mat>(kapa_pval[0]), join_rows(pval2, join_rows(conv_to<mat>::from(vars), rep(vars.n_elem,1,sela(dm)))+1));
        ntests+=vars.n_elem;
        pval(vars)=pmax<colvec>(pval(vars),pval2);
        ide = find(pval(vars) < alpha);
        vars = vars(ide);
        if ( vars.n_elem > 0  &&  max_k >= 2 ) {
            uvec sort_sela=sort(sela);
            for (int i=1;i<max_k;++i) {   //  Mixail comment:  change into this -> for ( i in 2:min( max_k, length(sela) ) ) {
                cand = find_combn<umat,uvec>(sort_sela, i+1); // anagkastika i+1 gia na vgalei 2 arithmous
                cand = cand.cols(find_arr_indices_c(cand==sela(dm)));  // Michail comment: delete this row
                for(int j=0;vars.n_elem > 0  &&  j < (int)cand.n_cols;++j){
                    devi_0 = calc_devi_0(y,Cbind(oness,x.cols(cand.col(j))),ini,maxiters,tol,test,lgmy,my,sy,d1,ylogy);
                    pval2 = add_term_c(y, Cbind(oness,x.cols(cand.col(j))),x.cols(vars),devi_0,ini,tol,true,parallel,maxiters,0);
                    //errors(482);
                    mat tmp1=conv_to<mat>::from(cand.col(j));
                    mat tmp=join_rows( pval2, join_rows(conv_to<mat>::from(vars), (mat)repmat(conv_to<mat>::from(cand.col(j).t()), vars.n_elem,1))+1);
                    //errors("row,col = ",tmp.n_rows,".",tmp.n_cols);
                    //errors("row,col = ",tmp1.n_rows,".",tmp1.n_cols);
                    kapa_pval[i] = join_cols( as<mat>(kapa_pval[i]), join_rows( pval2, join_rows(conv_to<mat>::from(vars), (mat)repmat(conv_to<rowvec>::from(cand.col(j)), vars.n_elem,1))+1));
                    ntests+=vars.n_elem;
                    pval(vars)=pmax<colvec>(pval(vars),pval2);
                    ide = find(pval(vars) < alpha);
                    vars = vars(ide);
                }  //## end  for ( j in 1:dim(cand)[2] ) {
            }  //## end  for ( i in 2:max_k ) {
        }  //## end  if ( max_k >= 2 ) {
        //errors("iiiiiiiiiiiiiiiiiii: ",vars.n_elem);
        colvec pval_vars=pval(vars);
        if(pval_vars.n_elem>0){
            sel=pval_vars.index_min();
            push_back(sela,vars(sel));
            vars = std_setdiff(vars, vars(sel));
        }
        //errors(505);
    } //## end  while ( vars.n_elem > 0  &  max_k >= 1 ) {

    /*if ( backward  && sela.n_elem > 0  ) {
    pv = pval(sela);
    sela = sela( order(pv) );
    bc = mmpcbackphase(y, x[, sela, drop = true], test,wei, max_k = max_k, threshold = threshold)
    met = bc$met
    sela = sela[met]
    pval[sela] = bc$pvalues
    ntests = ntests + bc$counter
    runtime = runtime + proc.time() - tic
} */
    timer.Stop();
    return List::create(_["selectedVars"] = sela+1, _["pvalues"] = pval, _["univ"] = univ,
                        _["kapa_pval"] = kapa_pval, _["max_k"] = max_k, _["threshold"] = alpha,
                        _["n.tests"] = ntests,_["runtime"]=timer.getTime());
}


RcppExport SEXP Rfast2_mmpc2(SEXP ySEXP,SEXP xSEXP,SEXP max_kSEXP,SEXP thresholdSEXP,SEXP testSEXP,SEXP Ini,SEXP parallelSEXP,SEXP maxitersSEXP,SEXP tolSEXP,SEXP backwardSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< vec >::type y(ySEXP);
    traits::input_parameter< mat >::type x(xSEXP);
    traits::input_parameter< int >::type max_k(max_kSEXP);
    traits::input_parameter< const double >::type threshold(thresholdSEXP);
    traits::input_parameter< const string >::type test(testSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const bool >::type backward(backwardSEXP);
    __result = mmpc2(y,x,max_k,threshold,test,Ini,parallel,maxiters,tol,backward);
    return __result;
END_RCPP
}
