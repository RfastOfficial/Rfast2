//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include <vector>
#include "Rfast2/templates.h"
#include "reg_lib_helper.h"
#include "mn.h"
#include "apply_funcs_templates.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

//[[Rcpp::plugins(openmp)]]
//[[Rcpp::plugins(cpp11)]]

void glm_poisson_2(mat x, vec y, const double lgmy,const double tol, const int maxiters, vec* ret){
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=pcols;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),m(n);
  mat L2,x_tr(n,pcols);
  double dif;
  b_old(0)=lgmy;
  x_tr=x.t();
  int i=2;
  for(dif=1.0;dif>tol;){
    yhat=x*b_old;
    m=(exp(yhat));
    L1=x_tr*(y-m);
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1,solve_opts::fast);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
    if(++i==maxiters)
      break;
  }

  ret[0] = m;
  ret[1] = b_old;
}

vec glm_logistic2(mat x, vec y,double *ini, vec yminmy, const double tol, const int maxiters){
  int n = x.n_rows, d = x.n_cols,i;

  vec est(n),m(n),p(n),w(n),be(d,fill::zeros);

  double d1 = ini[0], d2;

  be[0]= ini[1];

  mat der(d,1),der2(d,d),crossxy = cross_x_y<mat,mat,vec>(x, y);

  der = cross_x_y<mat,mat,vec>(x, yminmy);
  der2 = cross_x_y<mat,mat,vec>(x, x * ini[2]);
  be = be + solve(der2,der,solve_opts::fast);

  est = x*be;

  m = exp(-est);
  p = 1/(1+m);

  d2 = -2*calcDevRes(p,y,est);
  i=1;
  while(++i<maxiters && d1-d2>tol){

    d1 = d2;

    der = crossxy - cross_x_y<mat,mat,vec>(x, p);
    w = p%(1-p);

    der2 = cross_x_y<mat,mat,vec>(x, x.each_col()%w);

    be= be + solve(der2,der,solve_opts::fast);
    est = x*be;

    m = exp(-est);
    p = 1/(1+m);

    d2 = -2*calcDevRes(p,y,est);
  }

  return be;
}

double glm_logistic3(mat x, vec y,const double *ini, vec yminmy, const double tol, const int maxiters){
  int n = x.n_rows, d = x.n_cols,i;

  vec est(n),m(n),p(n),w(n),be(d,fill::zeros);

  double d1 = ini[0], d2;

  be[0]= ini[1];

  mat der(d,1),der2(d,d),crossxy = cross_x_y<mat,mat,vec>(x, y);

  der = cross_x_y<mat,mat,vec>(x, yminmy);
  der2 = cross_x_y<mat,mat,vec>(x, x * ini[2]);
  be = be + solve(der2,der,solve_opts::fast);

  est = x*be;

  m = exp(-est);
  p = 1/(1+m);

  d2 = -2*calcDevRes(p,y,est);
  i=1;
  while(++i<maxiters && d1-d2>tol){

    d1 = d2;

    der = crossxy - cross_x_y<mat,mat,vec>(x, p);
    w = p%(1-p);

    der2 = cross_x_y<mat,mat,vec>(x, x.each_col()%w);

    be= be + solve(der2,der,solve_opts::fast);
    est = x*be;

    m = exp(-est);
    p = 1/(1+m);

    d2 = -2*calcDevRes(p,y,est);
  }

  return d2;
}

vec logistic_only2(mat x, vec y, double *ini,vec yminmy,const double tol, const int maxiters){
  const unsigned int n=x.n_rows,pcols=x.n_cols;

  unsigned int i;
  int j;
  colvec be(2),expyhat(n),W(n,fill::zeros),x_col(n),x2_col(n),de(n);
  mat yhat(n,1),p(n,1);
  vec F(pcols);
  double d1,d2,t,dera2=0.0,sp=0.0,derb=0.0,dera=0.0,derab=0.0,derb2=0.0;
  const double sy=ini[3],lgmy=ini[1],d0 = ini[0];

  double W0=ini[2];
  double dera20=n*W0;
  colvec de0(n);
  de0=yminmy;
  double dera0=0;

  for(i=0;i<pcols;++i){
    d1=d0;
    be[0]=lgmy;
    be[1]=0;
    x_col=x.col(i);
    x2_col=arma::square(x_col);
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
    j=1;
    while(++j<maxiters && d1-d2>tol){
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

    F[i] = d2;
  }
  return F;
}

double glm_poisson3(mat x, vec y, const double lgmy, const double tol, int maxiters = 100){
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=pcols;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),m(n);
  mat L2,x_tr(n,pcols);
  double dif;
  b_old(0)=lgmy;
  x_tr=x.t();
  int i=2;
  for(dif=1.0;dif>tol;){
    yhat=x*b_old;
    m=(exp(yhat));
    L1=x_tr*(y-m);
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1,solve_opts::fast);
    dif=sum_with<std::abs,vec>(b_new-b_old);
    b_old=b_new;
    if(++i==maxiters)
      break;
  }

  //return 2.0*(ylogy-sum(y%yhat));
  return -2.0*sum(y%yhat);
}

vec weibull_mle2(vec x, int n, const double tol, const int maxiters){
  int i=2;
  vec lx = log(x),lx2 = lx%lx, y = x;
  double mlx = sum(lx)/n, co = sum(y%lx),sy = sum(y), fb = 1+mlx-co/sy,fb2 = -1 -(sum(y%lx2)*sy-co*co)/(sy*sy);
  double b1 = 1,b2 = 1 - fb/fb2;

  while (++i<maxiters && sum(abs(b2 - b1)) > tol) {
    b1 = b2;
    my_pow2(x,&y[0],b1,n);
    co = sum(y % lx);
    sy = sum(y);
    fb = 1/b1 + mlx - co/sy;
    fb2 = -1/(b1*b1) - (sum(y % lx2) * sy - co*co)/(sy*sy);
    b2 = b1 - fb/fb2;
  }

  vec param(3);
  double theta = pow(sy/n,1/b2);
  my_pow2(conv_to<vec>::from(x/theta),&y[0],b2,n);
  param[0] = b2;
  param[1] = theta;
  param[2] = n * log(b2) - n * b2 * log(theta) + (b2 - 1) * n * mlx - sum(y);

  return param;
}

double weib_reg2(vec y, mat x, vec ini, const double sly, const double tol, const int maxiters){
  int n = y.size(), d = x.n_cols;

  vec be(d,fill::zeros),com(n),lam;

  double ek = ini[0], lik1 = ini[2], lik2;
  be[0] = log(ini[1]);
  double yhat0 = exp(-be[0]);
  my_pow2(y*yhat0,&com[0],ek,n);

  rowvec sx = sum(x);

  vec logcom = log(com);
  vec comlogcom = com%logcom;

  double derk = n + ek * sly + ek*n*(-be[0]) - sum(comlogcom);
  double derk2 = derk - n - sum(comlogcom%logcom), k;
  k = log(ek) - derk/derk2;


  mat xcom = x.each_col()%com;
  rowvec derb =  sum(xcom)-sx;

  mat derb2 = -ek * cross_x_y<mat,mat,vec>(xcom, x);

  be = be -solve(derb2, derb.t(),solve_opts::fast);

  lam = -(x*be);

  vec yhat = exp(lam);
  ek = exp(k);
  my_pow2(y%yhat,&com[0],ek,n);

  lik2 = n*k+(ek-1)*sly+ek*sum(lam)-sum(com);
  int i = 1;

  while (++i<maxiters && lik2-lik1 > tol ) {
    lik1 = lik2;


    logcom = log(com);
    comlogcom = com%logcom;
    derk = n + ek * (sly + sum(lam)) - sum(comlogcom);
    derk2 = derk - n - sum(comlogcom%logcom);
    xcom = x.each_col()%com;
    derb =  sum(xcom)- sx;
    derb2 = -ek * cross_x_y<mat,mat,vec>(xcom, x);
    k = k - derk/derk2;
    be = be - solve(derb2, derb.t(),solve_opts::fast);
    lam = -(x*be);
    yhat = exp(lam);
    ek = exp(k);
    my_pow2(y%yhat,&com[0],ek,n);
    lik2 = n*k+(ek-1)*sly+ek*sum(lam)-sum(com);
  }

  return  lik2;
}

vec qpois_reg2(mat x, vec y,const double lgmy, const double ylogy,const double tol, const int maxiters){
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=pcols;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),m(n),phi(n);
  mat L2,x_tr(n,pcols);
  double dif;
  b_old(0)=lgmy;
  x_tr=x.t();
  int ij =2;
  for(dif=1.0;dif>tol;){
    yhat=x*b_old;
    m=(exp(yhat));
    phi=y-m;
    L1=x_tr*phi;
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1,solve_opts::fast);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
    if(++ij==maxiters)
      break;
  }
  double pheta = sum(arma::square(phi)/m)/(n-pcols);
  vec ret(3);
  ret(0) = 2.0*(ylogy-sum(y%yhat));
  ret(1) = pheta;
  ret(2) = (b_new[d-1]*b_new[d-1])/(pheta*((mat)arma::inv(L2))(d-1,d-1));
  return ret;
}

vec prop_reg2(mat x, vec y, const double *ini,vec yminmy, const double tol,const int maxiters){
  //ini = [d0,be0, my1minmy, ylogsums]
  const int n=x.n_rows,pcols=x.n_cols,d=pcols;
  colvec be(d,fill::zeros),yhat(n),expyhat,W(n,fill::zeros),p(n),u;

  double d1=ini[0],d2;

  be[0] = ini[1];

  mat der=cross_x_y<mat,mat,colvec>(x,yminmy);

  mat der2=cross_x_y<mat,mat,colvec>(x,x*ini[2]);
  be=be+solve(der2,der);
  yhat = x*be;
  expyhat=exp(-yhat);
  p = 1 / (1 + expyhat);
  d2=calcDevRes(p,y,expyhat);
  int i=2;
  for(;d2-d1>tol && i<maxiters;++i){
    d1=d2;
    u = y-p;
    der=cross_x_y<mat,mat,colvec>(x,u);
    W=p%(1-p);
    der2=cross_x_y<mat,mat,colvec>(x,x.each_col()%W);
    be=be+solve(der2,der);
    yhat = x*be;
    expyhat=exp(-yhat);
    p = 1 / (1 + expyhat);
    d2=calcDevRes(p,y,expyhat);
  }

  double dof = std::abs(n - d);
  double phi = sum((u % u)/(p%(1 - p)))/dof;
  mat slv = inv(der2);
  double  vb = phi * slv(d-1,d-1);

  vec ret(3);

  ret[0] = 2*(d2-ini[3]);
  ret[1] = phi;
  ret[2] = be[d-1]*be[d-1]/vb;

  return ret;
}

vec prop_regs2(mat x, vec y,double *ini, vec yminmy, const double tol, const int maxiters){
  const unsigned int n=x.n_rows,pcols=x.n_cols;

  unsigned int i=0;
  int ij=0;
  colvec be(2),expyhat(n),W(n,fill::zeros),x_col(n),x2_col(n),de(n), yhat(n), p(n);

  vec F(pcols);
  double d1=0,d2=0,t=0,dera2=0.0,sp=0.0,derb=0.0,dera=0.0,derab=0.0,derb2=0.0;
  const double lgmy=ini[1],d0 = ini[0],sy = ini[4];//d0 = -2*(sy*log(my)+(n-sy)*log(1-my));

  double W0=ini[2];
  double dera20=n*W0,vb;
  colvec de0(n);
  de0=yminmy;
  double dera0=0;

  for(i=0;i<pcols;++i){
    d1=d0;
    be[0]=lgmy;
    be[1]=0;
    x_col=x.col(i);
    x2_col=arma::square(x_col);
    derb=sum(de0%x_col);
    derab=W0*sum(x_col);
    derb2=W0*sum(x2_col);
    t=dera20 * derb2 - derab*derab;
    be[0]=be[0]+(derb2 * dera0 - derab * derb)/t;
    be[1]=be[1]+( - derab * dera0 + dera20 * derb )/t;

    yhat = be[0]+be[1]*x_col;
    expyhat=exp(-yhat);
    p = 1 / ( 1 + expyhat );

    d2 = calcDevRes(p,y,yhat);
    ij=1;

    while(++ij<maxiters && (d2-d1)>tol){
      d1=d2;
      W=p%(1-p);
      dera2=sum(W);
      sp=sum(p);

      de=y-p;
      dera=sy-sp;
      derb=sum(de%x_col);
      derab=sum(W%x_col);
      derb2=sum(W%x2_col);
      t=dera2 * derb2 - derab*derab;
      be[0]=be[0]+(derb2 * dera - derab * derb)/t;
      be[1]=be[1]+( - derab * dera + dera2 * derb )/t;
      yhat = be[0]+be[1]*x_col;
      expyhat=exp(-yhat);
      p = 1 / ( 1 + expyhat );

      d2 = calcDevRes(p,y,yhat);
    }

    vb=sum((arma::square(de))/W)/(n-2);
    vb=vb*dera2/(dera2*derb2-derab*derab);

    F[i]=-be[1]*be[1]/vb;
  }
  return F;
}

vec poisson_only2(mat x, vec y, int maxiters = 100){
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=2;
  unsigned int i;
  int ij;

  colvec b_old(d),b_new(d),L1(d),yhat(n);
  mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones);
  vec m(n),z_col_1(n);
  vec F(pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(mean(y));
  for(i=0;i<pcols;++i){
    b_old(0)=lgmeany;
    b_old(1)=0;
    z_col_1=x.col(i);
    z.col(1)=z_col_1;
    z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
    ij=2;
    for(dif=1.0;dif>0.000000001;){
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
    F[i]= -2.0*sum(y%yhat);
  }
  return F;
}

vec gold_rat3(double n, vec ni, vec ni2, double S, vec hi2,const int size, const double tol=1e-07){
  double a = 0, b = 50;
  const double ratio=0.618033988749895;
  double x1=b-ratio*b, x2=ratio*b;
  vec nix1 = ni*x1, nix2 = ni*x2, ni2hi2 = ni2%hi2;

  double f1 = calc_f(nix1, n, ni2hi2, S, x1, size);
  double f2 = calc_f(nix2, n, ni2hi2, S, x2, size);
  double bmina = b - a;
  while (abs(bmina)>tol){
    if(f2>f1){
      b=x2;
      bmina = b - a;
      x2=x1;
      f2=f1;
      x1=b - ratio * (bmina);
      nix1 = ni*x1;
      f1 = calc_f(nix1, n, ni2hi2, S, x1, size);
    }
    else {
      a=x1;
      bmina = b - a;
      x1=x2;
      f1=f2;
      x2=a + ratio * (bmina);
      nix2 = ni*x2;
      f2 = calc_f(nix2, n, ni2hi2, S, x2, size);
    }
  }
  vec ret(2);
  ret(0) = 0.5*(x1+x2);
  ret(1) = (f1+f2)/2;

  return ret;
}

vec rint_regs2(mat x, vec y, vec ni, IntegerVector id, int idmx, int idmn, vec sy, const double tol, const bool parallel, const int maxiters){
  int n = x.n_rows, D = x.n_cols;

  vec ni2 = ni%ni;
  vec my = sy/ni;
  double Sy = sum(sy);
  vec r = conv_to<vec>::from(cov(y,x));
  double mesi = Sy/n;
  vec xs = conv_to<vec>::from(sum(x));
  vec xs2 = conv_to<vec>::from(sum(x%x));
  vec vx = (xs2 - (xs%xs)/n)/(n - 1);
  vec b(D);
  b = r/vx;
  vec a(D);
  a = mesi - b % xs/n;
  mat be(D,2);
  be.col(0) = a;
  be.col(1) = b;
  vec stat(D);
  if(parallel){
#ifdef _OPENMP
#pragma omp parallel
{
#endif
  vec Xi(n), sxy(2), b1(2), tmpvec(n), tmpvec2(idmx), hi2(idmx),b2(2),B(2), mx(idmx);
  mat sx(idmx,2), temptcom(idmx,2), tcom(2,idmx), A(2,2);
  vec oneplnid;
  sx.col(0) = ni;
  sx.col(1) = ni;
  int ij=0;
  sxy[0] = Sy;
  double S=0,down=0,seb,se;
  vec d(2);
  mat  xx(2,2);
  xx(0,0) = n;
#ifdef _OPENMP
#pragma omp for
#endif
  for(int i = 0; i < D; ++i) {
    Xi = x.col(i);
    xx(0,1) = xs[i];
    xx(1,0) = xs[i];
    xx(1,1) = xs2[i];
    sx.col(1) = group_sum_helper<vec,vec,IntegerVector>(Xi, id, &idmn,&idmx);
    sxy[1] = sum(Xi % y);
    mx = sx.col(1)/ni;
    b1[0] = be.row(i)[0];
    b1[1] = be.row(i)[1];
    tmpvec = y - b1(0) - b1(1) * Xi;
    S = sum_with<square2<double>, vec>(tmpvec);
    tmpvec2 = my - b1(0) - b1(1) * mx;
    hi2 = tmpvec2%tmpvec2;
    d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
    oneplnid = (1+ ni * d[0]);
    temptcom.col(0) = sx.col(0)/oneplnid;
    temptcom.col(1) = sx.col(1)/oneplnid;
    tcom = -d[0] * temptcom.t();

    A = xx + tcom * sx;

    B = sxy + tcom * sy;

    down = A(0,0) * A(1,1) - A(0,1)*A(0,1);
    b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
    b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
    ij = 1;

    while(++ij<maxiters && std::abs(b1[0]+b1[1]-b2[0]-b2[1])>tol){
      b1 = b2;
      tmpvec = y - b1(0) - b1(1) * Xi;
      S = sum_with< square2<double>, vec>(tmpvec);
      tmpvec2 = my - b1(0) - b1(1) * mx;
      hi2 = tmpvec2%tmpvec2;
      d = gold_rat3(n, ni, ni2, S, hi2, idmx,tol);
      oneplnid = (1+ ni * d[0]);
      temptcom.col(0) = sx.col(0)/oneplnid;
      temptcom.col(1) = sx.col(1)/oneplnid;
      tcom = -d[0] * temptcom.t();
      A = xx + tcom * sx;
      B = sxy + tcom * sy;
      down = A(0,0) * A(1,1) - A(0,1) * A(0,1);
      b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
      b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
    }
    se = (S - d[0] * sum(ni2 % hi2/ oneplnid ) )/n;
    seb = A(0,0) / down * se;
    stat(i) = -b2(1)*b2(1)/ seb;
  }
#ifdef _OPENMP
}
#endif
  }
  else{
    vec Xi(n), sxy(2), b1(2), tmpvec(n), tmpvec2(idmx), hi2(idmx),b2(2),B(2), mx(idmx);
    mat sx(idmx,2), temptcom(idmx,2), tcom(2,idmx), A(2,2);
    vec oneplnid;
    sx.col(0) = ni;
    sx.col(1) = ni;
    int ij=0;
    sxy[0] = Sy;

    double S=0,down=0, seb,se;
    vec d(2);
    mat  xx(2,2);
    xx(0,0) = n;

    for(int i = 0; i < D; ++i) {
      Xi = x.col(i);
      xx(0,1) = xs[i];
      xx(1,0) = xs[i];
      xx(1,1) = xs2[i];
      sx.col(1) =  group_sum_helper<vec,vec,IntegerVector>(Xi, id, &idmn,&idmx);
      sxy[1] = sum(Xi % y);
      mx = sx.col(1)/ni;
      b1[0] = be.row(i)[0];
      b1[1] = be.row(i)[1];
      tmpvec = y - b1(0) - b1(1) * Xi;
      S = sum_with< square2<double>, vec>(tmpvec);
      //S = sum(tmpvec%tmpvec);
      tmpvec2 = my - b1(0) - b1(1) * mx;
      hi2 = tmpvec2%tmpvec2;
      d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
      oneplnid = (1+ ni * d[0]);
      temptcom.col(0) = sx.col(0)/oneplnid;
      temptcom.col(1) = sx.col(1)/oneplnid;
      tcom = -d[0] * temptcom.t();

      A = xx + tcom * sx;

      B = sxy + tcom * sy;

      down = A(0,0) * A(1,1) - A(0,1)*A(0,1);
      b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
      b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
      ij = 1;
      while(++ij<maxiters && std::abs(b1[0]+b1[1]-b2[0]-b2[1])>tol){
        //while(ij++<maxiters && sum(abs(b1-b2))>tol){
        b1 = b2;
        tmpvec = y - b1(0) - b1(1) * Xi;
        S = sum_with< square2<double>, vec>(tmpvec);
        //S = sum(tmpvec%tmpvec);
        tmpvec2 = my - b1(0) - b1(1) * mx;
        hi2 = tmpvec2%tmpvec2;
        d = gold_rat3(n, ni, ni2, S, hi2, idmx,tol);
        oneplnid = (1+ ni * d[0]);
        temptcom.col(0) = sx.col(0)/oneplnid;
        temptcom.col(1) = sx.col(1)/oneplnid;
        tcom = -d[0] * temptcom.t();
        A = xx + tcom * sx;
        B = sxy + tcom * sy;
        down = A(0,0) * A(1,1) - A(0,1) * A(0,1);
        b2(0) = (A(1,1) * B(0) - A(0,1) * B(1))/down;
        b2(1) = (- A(0,1) * B(0) + A(0,0) * B(1))/down;
      }
      se = (S - d[0] * sum(ni2 % hi2/ oneplnid ) )/n;
      seb = A(0,0) / down * se;
      stat(i) = -b2(1)*b2(1)/ seb;
    }
  }

  return stat;
}

double rint_reg2(mat x, vec y, vec ni, mat sx,vec sy, int idmx, const double tol, const int maxiters){
  int n = x.n_rows, p = x.n_cols;
  mat xx(p,p),sxy(p,1),mx;
  vec my;

  xx = cross_x<mat,mat>(x);

  sxy = cross_x_y<mat,mat,vec>(x,y);

  mx = sx.each_col()/ni;
  my = sy/ni;
  vec b1 = solve(xx,sxy,solve_opts::fast);
  vec tmp = y - x*b1;
  double S = sum(tmp%tmp);
  vec tmp2 = my-mx*b1;
  vec hi2 = tmp2%tmp2;
  vec ni2 = ni%ni;
  vec d(2);
  d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
  vec oneplnid = 1+ni*d(0);
  vec b2 = solve(xx - d(0)* cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid, sx), sxy -
    d(0) * cross_x_y<mat,mat,vec>(sx, sy/oneplnid),solve_opts::fast);
  int i = 2;

  while(++i<maxiters && apply_funcs<std::abs,mdiff<double>,double *,double *>(&b2[0],&b1[0],p,0) > tol) {
    b1 = b2;

    tmp = y - x*b1;
    S = sum_with<square2<double>,vec>(tmp);
    tmp2 = my-mx*b1;
    hi2 = tmp2%tmp2;

    d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
    oneplnid = 1+ni*d(0);
    b2 = solve(xx - d(0) * cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid, sx), sxy -
      d(0) * cross_x_y<mat,mat,vec>(sx, sy/oneplnid),solve_opts::fast);
  }

  double info2 = (S-d(0)*sum(ni2%hi2/oneplnid))/n;

  double belast = b2(p-1);
  vec seb = ((mat)arma::inv(xx-d(0)*cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid,sx))).diag()*info2;

  return belast*belast/seb(p-1);
}

double spml_reg2(mat pu, mat x, const double tol, const int maxiters){
  mat tx = x.t(), XX = solve(tx*x,tx,solve_opts::fast);
  int n = x.n_rows, D = x.n_cols;
  mat u = mat(pu.begin(),n,2,false), be = XX*u, mu = x*be;

  vec tau = sum(u%mu,1), ptau = pnormc(tau), ci2(pu.begin_col(2),n,false), cisi(pu.begin_col(3),n,false), si2(pu.begin_col(4),n,false);

  double f= -0.5, con = 2.506628274631, lik1 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n), lik2;
  vec rat, psit2, psit = tau + ptau / ( exp(f * tau%tau)/con + tau % ptau );
  be = XX * (u.each_col() % psit);
  mu = mu.each_col()% psit;
  tau = sum(u%mu,1);
  ptau = pnormc(tau);
  lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);
  mat der2, a11,a12,a22;

  mat der, slv;
  int i = 2;

  while(++i<maxiters && lik2 - lik1 > tol){
    lik1 = lik2;

    rat = ptau / ( exp(f * tau%tau)/con + tau % ptau );
    psit = tau + rat;

    psit2 = 2 -((tau+rat)%rat);
    der = cross_x_y<mat,mat,vec>(x, u.each_col()%psit-mu);

    der.reshape(2*D,1);



    a11 = cross_x_y<mat,mat,vec>(x, x.each_col() % (psit2%ci2 - 1));
    a12 = cross_x_y<mat,mat,vec>(x, x.each_col() % (psit2%cisi));
    a22 = cross_x_y<mat,mat,vec>(x, x.each_col() % (psit2%si2 - 1));

    der2 = join_cols(join_rows(a11,a12),join_rows(a12,a22));

    slv = solve(der2, der,solve_opts::fast);
    slv.reshape(D,2);
    be = be - slv;

    mu = x * be;
    tau = sum(u%mu,1);
    ptau = pnormc(tau);

    lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);
  }

  return lik2 - n*1.8378770664;
}


double multinom_reg2(mat y, mat x, mat e, rowvec m0, rowvec b0, const double tol, const int maxiters){
  int n = x.n_rows, p = x.n_cols, d = y.n_cols,  pd = p*d;

  mat b1(p,d,fill::zeros);

  b1.row(0) = b0;

  rowvec der(pd);

  mat der2(pd,pd,fill::zeros), crossress, b2(p,d), m1, m;

  double exp20 = exp(20);

  vec slv;
  mat xCrossX = cross_x_y<mat,mat,vec>(x,x);

  int pi,j,pj, ij;
  arma::span spani, spanj;
  mat::iterator slvit, b2it = b2.begin(), b1it;
  for(int i = 0; i<d; ++i){
    pi = i*p;
    spani = span(pi,p-1+pi);

    der(spani) = sum(x.each_col()%e.col(i));

    for(j = i; j < d; ++j){
      pj = j*p;
      spanj = span(pj,p-1+pj);
      if(i!=j){
        crossress = -(m0[i] * m0[j]) * xCrossX;
        der2(spani, spanj) = crossress;
        der2(spanj, spani) = crossress;
      }
      else{
        crossress = (m0[i] * (1 - m0[i])) * xCrossX;
        der2(spani, spani) = crossress;
      }
    }
  }

  slv = solve(der2, der.t(),solve_opts::fast);

  b1it = b1.begin();
  slvit = slv.begin();

  apply_funcs<madd<double>, mat::iterator, vec::iterator, mat::iterator>(b1it,slvit, b2it, pd);

  ij=1;
  colvec one(n,fill::ones);

  while(++ij<maxiters && apply_funcs<abs, mdiff<double>,
        mat::iterator, mat::iterator>(b1it,b2it,pd,0)> tol) {
    b1 = b2;
    m1 = clamp(exp(x*b1),0,exp20);

    m = m1.each_col()/ (sum(m1,1) + 1);

    e = y - m;

    for(int i = 0; i<d; ++i){
      pi = i*p;
      spani = span(pi,p-1+pi);
      der(spani) = sum(x.each_col()%e.col(i));

      for(j = i; j < d; ++j){
        pj = j*p;
        spanj = span(pj,p-1+pj);
        if(i!=j){
          crossress = -cross_x_y<mat,mat,vec>(x.each_col() % (m.col(i) % m.col(j)), x);
          der2(spani, spanj) = crossress;
          der2(spanj, spani) = crossress;
        }
        else{
          crossress = cross_x_y<mat,mat,vec>(x.each_col() % (m.col(i) % (one-m.col(i))), x);
          der2(spani, spani) = crossress;
        }
      }
    }

    slv = solve(der2, der.t(),solve_opts::fast);

    b1it = b1.begin();
    slvit = slv.begin();

    apply_funcs<madd<double>, mat::iterator, vec::iterator, mat::iterator>(b1it,slvit, b2it, pd);
  }

  mat m2 = join_horiz(one,m1);
  m2 = m2.each_col()/sum(m2,1);

  return -mreg_loglic(y,m2);
}

vec normlog_reg2(vec y,mat x, vec lgy01, const double tol, const int maxiters){
  int n = x.n_rows, d = x.n_cols, elems = n*d;
  mat xt = x.t();
  colvec b1 = solve(xt * x, xt * lgy01,solve_opts::fast);
  colvec yhat = (exp(x * b1));

  mat tmpmat(n,d);


  mat com = x.each_col()%(yhat % yhat);

  mat com2 = x.each_col()%(yhat % y);

  apply_funcs<mdiff<double>, double *,double *, double *>(&com[0],&com2[0],&tmpmat[0],elems);


  vec der = sum(tmpmat).t();

  apply_funcs<madd<double>, double *,double *, double *>(&tmpmat[0],&com[0],&tmpmat[0],elems);
  mat der2 = xt*tmpmat;

  vec b2 = b1 - solve(der2, der,solve_opts::fast);
  int i = 1;

  while(++i<maxiters && apply_funcs<std::abs,mdiff<double>,double *,double *>(&b1[0],&b2[0],d,0) > tol){
    b1 = b2;
    yhat = (exp(x * b1));
    com = x.each_col()%(yhat % yhat);
    com2 = x.each_col()%(yhat % y);
    apply_funcs<mdiff<double>, double *,double *, double *>(&com[0],&com2[0],&tmpmat[0],elems);

    der = sum(tmpmat).t();

    apply_funcs<madd<double>, double *,double *, double *>(&tmpmat[0],&com[0],&tmpmat[0],elems);
    der2 = xt*tmpmat;

    b2 = b1 - solve(der2, der,solve_opts::fast);
  }
  double deviance = apply_funcs<square2<double>, mdiff<double>, double *, double *>(&y[0], &yhat[0], n, 0);

  vec ret(2);
  ret(0) =  deviance;
  ret(1) = deviance/(n-d);
  return ret;
}

double vmf_mle2(double nR, const int n, const double tol, const double maxiters){
  double apk, k2, r = nR/n, R2 = r*r, k1 = r * (2 - R2)/(1 - R2);

  int i = 1;
  if(k1 < 1e+05){
    apk = R::bessel_i(k1,1, 1)/R::bessel_i(k1,0, 1);
    k2 = k1 - (apk - r)/(1 - apk*apk - 1/k1 * apk);
    while (++i<maxiters && std::abs(k2 - k1) > tol) {
      k1 = k2;
      apk = R::bessel_i(k1,1,1)/R::bessel_i(k1,0,1);
      k2 = k1 - (apk - r)/(1 - apk*apk - 1/k1 * apk);
    }
    return k2;
  }

  return k1;
}

double spml_mle2(mat u, vec ci2, vec cisi, vec si2, const int n, const double tol, const int maxiters){
  vec su(2);
  su(0) = sum(u.col(0)), su(1) = sum(u.col(1));

  double nR = sqrt(sum(su%su)), kappa = vmf_mle2(nR, n, tol, maxiters);

  vec mu = su*(1/nR);

  vec mu1 = mu*kappa;
  double f = -0.5, con = 2.506628274631;
  vec tau = u*mu1, ptau = pnormc(tau);

  vec rat = ptau/(exp(f * tau%tau)/con + tau % ptau);

  vec psit = tau + rat;
  vec psit2 = 2 - rat%(tau + rat);

  vec der(2);
  der[0] = sum(u.col(0)%psit) - n * mu1[0];
  der[1] = sum(u.col(1)%psit) - n * mu1[1];

  double dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;

  double down = dera2 * derb2 - derab*derab;
  vec mu2(2);
  mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
  mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

  int i = 1;
  while (++i<maxiters && sum(abs(mu2 - mu1)) > tol) {
    mu1 = mu2;
    tau = u*mu1;
    ptau = pnormc(tau);
    rat = ptau/(exp(f * (tau%tau))/con + tau % ptau);
    psit = tau + rat;
    psit2 = 2 - rat%(tau + rat);
    der[0] = sum(u.col(0)%psit) - n * mu1[0];
    der[1] = sum(u.col(1)%psit) - n * mu1[1];
    dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;
    down = dera2 * derb2 - derab*derab;
    mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
    mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

  }

  return -0.5 * n * (mu2[0]*mu2[0]+mu2[1]*mu2[1]) + sum_with<log1p, colvec>(((tau % ptau) * con)/exp(f*tau%tau)) - n * 1.83787706640935;
}

vec spml_regs2(mat pu, mat x, const double tol, const int maxiters, const bool parallel){
  int n = x.n_rows, D = x.n_cols;

  vec ci2(pu.begin_col(2),n,false), cisi(pu.begin_col(3),n,false), si2(pu.begin_col(4),n,false);
  mat u(pu.begin(),n,2,false);
  double f = -0.5, con = 2.506628274631;

  vec one(n,fill::ones);
  double ini = spml_mle2(u,ci2,cisi,si2,n,tol,maxiters);

  vec ret(D);

  if(parallel){
#ifdef _OPENMP
#pragma omp parallel
{
#endif
  double sumX, lik1, lik2;
  int i;
  mat xx(2,2,fill::zeros),X(n,2), B,mu,a11,a12,a22,der2(4,4),XX,der(4,1);
  vec tau,ptau,psit,rat,psit2,slv;
  mat::iterator tmpit = der2.begin(),it;
  X.col(0) = one;

  xx(0) = n;
#ifdef _OPENMP
#pragma omp for
#endif
  for(int j = 0; j < D; ++j){
    X.col(1) = x.col(j);
    sumX = sum(X.col(1));
    xx(1) = sumX;
    xx(2) = sumX;
    xx(3) = sum(X.col(1)%X.col(1));
    XX = solve(xx, X.t());
    B =  XX * u;
    mu = X * B;
    tau = sum(u % mu,1);
    ptau = pnormc(tau);

    lik1 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);

    rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
    psit = tau + rat;
    psit2 = 2 -((tau+rat)%rat);
    a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
    der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
    a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
    a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
    a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
    it = tmpit;

    // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
    (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
    (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
    (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
    (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

    slv = solve(der2,der);
    B[0] = B[0]-slv[0];
    B[1] = B[1]-slv[1];
    B[2] = B[2]-slv[2];
    B[3] = B[3]-slv[3];

    mu = X*B;

    tau = sum(u % mu,1);
    ptau = pnormc(tau);

    lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);

    i=2;
    while(++i<maxiters && (lik2 - lik1) > tol){
      lik1 = lik2;
      rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
      psit = tau + rat;
      psit2 = 2 - tau % rat - rat%rat;
      a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
      der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
      a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
      a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
      a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
      it = tmpit;

      // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
      (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
      (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
      (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
      (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

      slv = solve(der2,der);
      B[0] = B[0]-slv[0];
      B[1] = B[1]-slv[1];
      B[2] = B[2]-slv[2];
      B[3] = B[3]-slv[3];

      mu = X*B;

      tau = sum(u % mu,1);
      ptau = pnormc(tau);

      lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);
    }

    ret(j) = 2 * (lik2 - n * 1.83787706640935 - ini);
  }
#ifdef _OPENMP
}
#endif
  }
  else{
    double sumX, lik1, lik2;
    int i;
    mat xx(2,2,fill::zeros),X(n,2), B,mu,a11,a12,a22,der2(4,4),XX,der(4,1);
    vec tau,ptau,psit,rat,psit2,slv;
    mat::iterator tmpit = der2.begin(),it;
    X.col(0) = one;

    xx(0) = n;

    for(int j = 0; j < D; ++j){
      X.col(1) = x.col(j);
      sumX = sum(X.col(1));
      xx(1) = sumX;
      xx(2) = sumX;
      xx(3) = sum(X.col(1)%X.col(1));
      XX = solve(xx, X.t());
      B =  XX * u;
      mu = X * B;
      tau = sum(u % mu,1);
      ptau = pnormc(tau);

      lik1 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);

      rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
      psit = tau + rat;
      psit2 = 2 - tau % rat - rat%rat;
      a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
      der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
      a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
      a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
      a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
      it = tmpit;

      // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
      (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
      (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
      (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
      (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

      slv = solve(der2,der);
      B[0] = B[0]-slv[0];
      B[1] = B[1]-slv[1];
      B[2] = B[2]-slv[2];
      B[3] = B[3]-slv[3];

      mu = X*B;

      tau = sum(u % mu,1);
      ptau = pnormc(tau);

      lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);

      i=2;
      while(++i<maxiters && (lik2 - lik1) > tol){
        lik1 = lik2;
        rat = ptau / ( exp(f * (tau%tau))/con + tau%ptau );
        psit = tau + rat;
        psit2 = 2 - tau % rat - rat%rat;
        a11 = cross_x_y<mat,mat,vec>(X, u.each_col()%psit-mu);
        der[0] = a11[0],der[1] = a11[1],der[2] = a11[2],der[3] = a11[3];
        a11 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%ci2 - 1));
        a12 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%cisi));
        a22 = cross_x_y<mat,mat,vec>(X, X.each_col() % (psit2%si2 - 1));
        it = tmpit;

        // in our case equivalent to cbind( rbind(a11, a12), rbind(a12, a22) )
        (*(it)) = a11[0], (*(++it)) = a11[2], (*(++it)) = a12[0], (*(++it)) = a12[2];
        (*(++it)) = a11[1], (*(++it)) = a11[3], (*(++it)) = a12[1], (*(++it)) = a12[3];
        (*(++it)) = a12[0], (*(++it)) = a12[2], (*(++it)) = a22[0], (*(++it)) = a22[2];
        (*(++it)) = a12[1], (*(++it)) = a12[3], (*(++it)) = a22[1], (*(++it)) = a22[3];

        slv = solve(der2,der);
        B[0] = B[0]-slv[0];
        B[1] = B[1]-slv[1];
        B[2] = B[2]-slv[2];
        B[3] = B[3]-slv[3];

        mu = X*B;

        tau = sum(u % mu,1);
        ptau = pnormc(tau);

        lik2 = calc_spml_loglik(mu.begin_col(0),mu.begin_col(1),&tau[0],&ptau[0],n);
      }
      ret(j) = 2 * (lik2 - n * 1.83787706640935 - ini);
    }
  }

  return ret;
}

vec multinom_regs2(mat Y, mat x, const double tol, const bool parallel, const int maxiters){
  int n = x.n_rows, D = x.n_cols;

  //mat Y1 = design_matrix_helper<mat,NumericVector>(Y0);

  vec poia = indexesOfNum(Y,1);
  int poiasize = poia.size();

  rowvec m0 = mean(Y);

  double ini = calc_multinom_ini(Y,conv_to<vec>::from(m0));

  Y.shed_col(0);

  n = Y.n_rows;
  int d = Y.n_cols;

  mat b10(2,d,fill::zeros),e0;

  rowvec b0(d);
  m0.shed_col(0);
  /*for(int i = 0; i < d; i++)
   m0[i] = tmpvec[i+1];*/

  b0 = log(m0/(1-m0));
  b10.row(0) = b0;
  e0 = Y.each_row()-m0;

  mat id = create_id_mat(d);

  double exp20 = exp(20);

  int dx2 = 2*d;
  vec ret(D);
  colvec one(n,fill::ones);
  /*--------------------------------------------------------*/
  if(parallel){
#ifdef _OPENMP
#pragma omp parallel
{
#endif
  mat dera(n,dx2),der2(dx2,dx2,fill::zeros), b1,b2(2,d),m1,m,e1,crossress,X(n,2), xCrossx;
  vec der(dx2),idcoli,idcolj,slv;
  mat::iterator slvit, b2it, b1it;
  int i=0,j=0,ij=0;
  X.col(0) = one;
#ifdef _OPENMP
#pragma omp for
#endif
  for(int l = 0; l < D; ++l) {
    X.col(1) = x.col(l);
    xCrossx = cross_x_y<mat,mat,vec>(X, X);
    for(i = 0; i < d; ++i) {
      idcoli = id.col(i);

      dera.col(idcoli[0]) = e0.col(i)%X.col(0);
      dera.col(idcoli[1]) = X.col(1)%e0.col(i);
      for (j = i; j < d; ++j) {
        if (i != j) {
          idcolj = id.col(j);
          crossress = -(m0(i) * m0(j))*xCrossx;
          der2(idcolj[0], idcoli[0]) =  crossress(0,0);
          der2(idcolj[0], idcoli[1]) =  crossress(0,1);
          der2(idcolj[1], idcoli[0]) =  crossress(1,0);
          der2(idcolj[1], idcoli[1]) =  crossress(1,1);

          der2(idcoli[0], idcolj[0]) =  crossress(0,0);
          der2(idcoli[0], idcolj[1]) =  crossress(0,1);
          der2(idcoli[1], idcolj[0]) =  crossress(1,0);
          der2(idcoli[1], idcolj[1]) =  crossress(1,1);
        }
        else {
          crossress = ((m0(i) * (1 - m0(i)))) * xCrossx;

          der2(idcoli[0], idcoli[0]) =  crossress(0,0);
          der2(idcoli[0], idcoli[1]) =  crossress(0,1);
          der2(idcoli[1], idcoli[0]) =  crossress(1,0);
          der2(idcoli[1], idcoli[1]) =  crossress(1,1);
        }
      }
    }
    der = conv_to<vec>::from(sum(dera));

    b1 = b10;

    slv = solve(der2, der);

    b2it = b2.begin();
    b1it = b1.begin();
    slvit = slv.begin();
    apply_funcs<madd<double>, mat::iterator, mat::iterator, mat::iterator>(b1it,slvit, b2it, dx2);

    ij=2;
    while(++ij<maxiters && apply_funcs<abs, mdiff<double>,
          mat::iterator, mat::iterator>(b1it,b2it,dx2,0) > tol) {
      b1 = b2;

      m1 = clamp(exp(X*b1),0,exp20);

      m = m1.each_col()/ (sum(m1,1) + 1);

      e1 = Y - m;
      for(i = 0; i<d; ++i) {
        idcoli = id.col(i);
        dera.col(idcoli[0]) = e1.col(i)%X.col(0);
        dera.col(idcoli[1]) = e1.col(i)%X.col(1);

        for (j = 0; j<d; ++j) {
          if (i != j) {
            idcolj = id.col(j);

            crossress = -cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % m.col(j)), X);
            der2(idcoli[0], idcolj[0]) =  crossress(0,0);
            der2(idcoli[0], idcolj[1]) =  crossress(0,1);
            der2(idcoli[1], idcolj[0]) =  crossress(1,0);
            der2(idcoli[1], idcolj[1]) =  crossress(1,1);

            der2(idcolj[0], idcoli[0]) =  crossress(0,0);
            der2(idcolj[0], idcoli[1]) =  crossress(0,1);
            der2(idcolj[1], idcoli[0]) =  crossress(1,0);
            der2(idcolj[1], idcoli[1]) =  crossress(1,1);
          }
          else {
            crossress = cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % (one - m.col(i))), X);

            der2(idcoli[0], idcoli[0]) =  crossress(0,0);
            der2(idcoli[0], idcoli[1]) =  crossress(0,1);
            der2(idcoli[1], idcoli[0]) =  crossress(1,0);
            der2(idcoli[1], idcoli[1]) =  crossress(1,1);
          }
        }
      }
      der = conv_to<vec>::from(sum(dera));

      slv = solve(der2, der);
      b2it = b2.begin();
      b1it = b1.begin();
      slvit = slv.begin();
      apply_funcs<madd<double>, mat::iterator, mat::iterator, mat::iterator>(b1it,slvit, b2it, dx2);
    }

    m1.insert_cols(0,one);
    m1 = m1.each_col() / sum(m1,1);

    ret(l) = 2 * calcSumLog(m1,poia,poiasize) - ini;
  }
#ifdef _OPENMP
}
#endif
  }
  else{
    mat dera(n,dx2),der2(dx2,dx2,fill::zeros), b1,b2(2,d),m1,m,e1,crossress,X(n,2),xCrossx;
    vec der(dx2),idcoli,idcolj,slv;
    mat::iterator slvit, b2it, b1it;
    int i=0,j=0,ij=0;
    X.col(0) = one;

    for(int l = 0; l < D; ++l) {
      X.col(1) = x.col(l);
      xCrossx = cross_x_y<mat,mat,vec>(X, X);
      for(i = 0; i < d; ++i) {
        idcoli = id.col(i);

        dera.col(idcoli[0]) = e0.col(i)%X.col(0);
        dera.col(idcoli[1]) = X.col(1)%e0.col(i);
        for (j = i; j < d; ++j) {
          if (i != j) {
            idcolj = id.col(j);
            crossress = -(m0(i) * m0(j))*xCrossx;
            der2(idcolj[0], idcoli[0]) =  crossress(0,0);
            der2(idcolj[0], idcoli[1]) =  crossress(0,1);
            der2(idcolj[1], idcoli[0]) =  crossress(1,0);
            der2(idcolj[1], idcoli[1]) =  crossress(1,1);

            der2(idcoli[0], idcolj[0]) =  crossress(0,0);
            der2(idcoli[0], idcolj[1]) =  crossress(0,1);
            der2(idcoli[1], idcolj[0]) =  crossress(1,0);
            der2(idcoli[1], idcolj[1]) =  crossress(1,1);
          }
          else {
            crossress = ((m0(i) * (1 - m0(i)))) * xCrossx;

            der2(idcoli[0], idcoli[0]) =  crossress(0,0);
            der2(idcoli[0], idcoli[1]) =  crossress(0,1);
            der2(idcoli[1], idcoli[0]) =  crossress(1,0);
            der2(idcoli[1], idcoli[1]) =  crossress(1,1);
          }
        }
      }
      der = conv_to<vec>::from(sum(dera));

      b1 = b10;

      slv = solve(der2, der);

      b2it = b2.begin();
      b1it = b1.begin();
      slvit = slv.begin();
      apply_funcs<madd<double>, mat::iterator, mat::iterator, mat::iterator>(b1it,slvit, b2it, dx2);

      ij=2;
      while(++ij<maxiters && apply_funcs<abs, mdiff<double>,
            mat::iterator, mat::iterator>(b1it,b2it,dx2,0) > tol) {
        b1 = b2;

        m1 = clamp(exp(X*b1),0,exp20);

        m = m1.each_col()/ (sum(m1,1) + 1);

        e1 = Y - m;
        for(i = 0; i<d; ++i) {
          idcoli = id.col(i);
          dera.col(idcoli[0]) = e1.col(i)%X.col(0);
          dera.col(idcoli[1]) = e1.col(i)%X.col(1);

          for (j = 0; j<d; ++j) {
            if (i != j) {
              idcolj = id.col(j);

              crossress = -cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % m.col(j)), X);
              der2(idcoli[0], idcolj[0]) =  crossress(0,0);
              der2(idcoli[0], idcolj[1]) =  crossress(0,1);
              der2(idcoli[1], idcolj[0]) =  crossress(1,0);
              der2(idcoli[1], idcolj[1]) =  crossress(1,1);

              der2(idcolj[0], idcoli[0]) =  crossress(0,0);
              der2(idcolj[0], idcoli[1]) =  crossress(0,1);
              der2(idcolj[1], idcoli[0]) =  crossress(1,0);
              der2(idcolj[1], idcoli[1]) =  crossress(1,1);
            }
            else {
              crossress = cross_x_y<mat,mat,vec>(X.each_col()%(m.col(i) % (one - m.col(i))), X);

              der2(idcoli[0], idcoli[0]) =  crossress(0,0);
              der2(idcoli[0], idcoli[1]) =  crossress(0,1);
              der2(idcoli[1], idcoli[0]) =  crossress(1,0);
              der2(idcoli[1], idcoli[1]) =  crossress(1,1);
            }
          }
        }
        der = conv_to<vec>::from(sum(dera));

        slv = solve(der2, der);
        b2it = b2.begin();
        b1it = b1.begin();
        slvit = slv.begin();
        apply_funcs<madd<double>, mat::iterator, mat::iterator, mat::iterator>(b1it,slvit, b2it, dx2);
      }

      m1.insert_cols(0,one);
      m1 = m1.each_col() / sum(m1,1);

      ret(l) = 2 * calcSumLog(m1,poia,poiasize) - ini;
    }
  }

  return ret;
}

vec weib_regs2(vec y, mat x, vec ini, const double sly, const double tol, const int maxiters, const bool parallel){
  int n = y.size(), d = x.n_cols;

  vec one(n,fill::ones);


  vec com0(n), logcom0, comlogcom0;

  double be0 = log(ini[1]), ek0 = ini[0], lik0 = ini[2], yhat0 = exp(-be0),k0;
  my_pow2(y*yhat0,&com0[0],ek0,n);

  logcom0 = log(com0);
  comlogcom0 = com0%logcom0;

  double derk0 = n + ek0 * sly + ek0*n*(-be0) - sum(comlogcom0);
  double derk20 = derk0 - n - sum(comlogcom0%logcom0);
  k0 = log(ek0) - derk0/derk20;

  vec ret(d);
  /*---------------------------------------*/
  if(parallel){
#ifdef _OPENMP
#pragma omp parallel
{
#endif
  vec derb(2),derb2(3),xcolicom,slv(2), be(2), com,lam,logcom,comlogcom,yhat;
  mat tmpX(n,2);//,xcom(n,2);
  double scom,sxcolicom, sxcoli, ek,lik1,lik2,k,derk,derk2;
  tmpX.col(0) = one;
  int iters;
#ifdef _OPENMP
#pragma omp for
#endif
  for(int i = 0; i < d; ++i){
    ek = ek0, lik1 = lik0, k=k0, derk = derk0, derk2 = derk20;
    com=com0,logcom = logcom0,comlogcom=comlogcom0;

    be[0] = be0;
    be[1] = 0;

    tmpX.col(1) = x.col(i);
    sxcoli = sum(tmpX.col(1));

    scom = sum(com);
    xcolicom = tmpX.col(1)%com;
    sxcolicom = sum(xcolicom);

    derb[0] = scom-n;
    derb[1] = sxcolicom - sxcoli;


    derb2(0) = -scom *ek;
    derb2(1) = -sum(xcolicom) * ek;
    derb2(2) = -sum(tmpX.col(1)%xcolicom) * ek;

    slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
    slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
    be = be - slv;

    lam = -(tmpX*be);
    yhat = exp(lam);
    ek = exp(k);

    my_pow2(y%yhat,&com[0],ek,n);
    scom = sum(com);
    lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;

    iters = 2;

    while (++iters<maxiters && lik2-lik1 > tol ) {
      lik1 = lik2;

      logcom = log(com);
      comlogcom = com%logcom;
      derk = n + ek * (sly + sum(lam)) - sum(comlogcom);
      derk2 = derk - n - sum(comlogcom%logcom);

      xcolicom = tmpX.col(1)%com;
      sxcolicom = sum(xcolicom);

      derb[0] = scom-n;
      derb[1] = sxcolicom - sxcoli;

      derb2[0] = -scom *ek;
      derb2[1] = -sum(xcolicom) * ek;
      derb2[2] = -sum(tmpX.col(1)%xcolicom) * ek;

      slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
      slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
      be = be - slv;
      k = k - derk/derk2;

      lam = -(tmpX*be);
      yhat = exp(lam);
      ek = exp(k);

      my_pow2(y%yhat,&com[0],ek,n);
      scom = sum(com);
      lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;
    }
    ret(i) = 2*(lik2-lik0);
  }
#ifdef _OPENMP
}
#endif
  }
  else{
    vec derb(2),derb2(3),xcolicom, be(2),slv(2), com,lam,logcom,comlogcom,yhat;
    mat tmpX(n,2);
    double scom,sxcolicom, sxcoli, ek,lik1,lik2,k,derk,derk2;
    int iters;
    tmpX.col(0) = one;
    for(int i = 0; i < d; ++i){
      ek = ek0, lik1 = lik0, k=k0, derk = derk0, derk2 = derk20;
      com=com0,logcom = logcom0,comlogcom=comlogcom0;

      be[0] = be0;
      be[1] = 0;

      tmpX.col(1) = x.col(i);
      sxcoli = sum(tmpX.col(1));

      scom = sum(com);
      xcolicom = tmpX.col(1)%com;
      sxcolicom = sum(xcolicom);

      derb[0] = scom-n;
      derb[1] = sxcolicom - sxcoli;


      derb2(0) = -scom *ek;
      derb2(1) = -sum(xcolicom) * ek;
      derb2(2) = -sum(tmpX.col(1)%xcolicom) * ek;

      slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
      slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
      be = be - slv;

      lam = -(tmpX*be);
      yhat = exp(lam);
      ek = exp(k);

      my_pow2(y%yhat,&com[0],ek,n);
      scom = sum(com);
      lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;

      iters = 2;

      while (++iters<maxiters && lik2-lik1 > tol ) {
        lik1 = lik2;

        logcom = log(com);
        comlogcom = com%logcom;
        derk = n + ek * (sly + sum(lam)) - sum(comlogcom);
        derk2 = derk - n - sum(comlogcom%logcom);

        xcolicom = tmpX.col(1)%com;
        sxcolicom = sum(xcolicom);

        derb[0] = scom-n;
        derb[1] = sxcolicom - sxcoli;

        derb2[0] = -scom *ek;
        derb2[1] = -sum(xcolicom) * ek;
        derb2[2] = -sum(tmpX.col(1)%xcolicom) * ek;

        slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
        slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
        be = be - slv;
        k = k - derk/derk2;

        lam = -(tmpX*be);
        yhat = exp(lam);
        ek = exp(k);

        my_pow2(y%yhat,&com[0],ek,n);
        scom = sum(com);
        lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;
      }
      ret(i) = 2*(lik2-lik0);
    }
  }



  return ret;
}

vec normlog_regs2(vec y,mat x, vec ly, const double tol,const bool parallel,const int maxiters){
  int n = x.n_rows;
  int D = x.n_cols;

  //ly = log(y+0.1);

  int con = var(y)*(n-1);

  double my = sum(ly)/n;

  vec ret(D);
  if(parallel) {
#ifdef _OPENMP
#pragma omp parallel
{
#endif
  vec tmpX, aold(2), tmpX2, yhat, com, com2, anew(2);
  double dera, dera2, derb, derb2, derab, com3, sumcom2, divider, va;
  rowvec berow(2);
  int ij;
#ifdef _OPENMP
#pragma omp for
#endif
  for(int j = 0; j<D; ++j){
    tmpX = x.col(j);

    aold[1] = (cov(ly,tmpX)/var(tmpX))[0];
    aold[0] = my - aold[1]*sum(tmpX)/n;

    tmpX2 = tmpX % tmpX;

    yhat = exp(aold[0] + aold[1] * tmpX);
    com = y % yhat;
    com2 = yhat % yhat;
    sumcom2 = sum(com2);
    com3 = sum(com2 % tmpX);
    dera = sum(com) - sumcom2;
    dera2 = dera - sumcom2;
    derb = sum(com % tmpX) - com3;
    derb2 = sum(com % tmpX2) - 2 * sum(com2 % tmpX2);
    derab = derb - com3;

    divider = dera2 * derb2 - derab *derab;

    anew[0] = aold[0] - (derb2 * dera - derab * derb)/divider;
    anew[1] = aold[1] - (-derab * dera + dera2 * derb)/divider;

    ij =2;
    while (++ij<maxiters && std::abs(anew[0]+anew[1] - aold[0]-aold[1])  > tol ){
      aold = anew;

      yhat = exp(aold[0] + aold[1] * tmpX);
      com = y % yhat;
      com2 = yhat % yhat;
      sumcom2 = sum(com2);
      com3 = sum(com2 % tmpX);
      dera = sum(com) - sumcom2;
      dera2 = dera - sumcom2;
      derb = sum(com % tmpX) - com3;
      derb2 = sum(com % tmpX2) - 2*sum(com2 % tmpX2);
      derab = derb - com3;

      divider = dera2 * derb2 - derab *derab;

      anew[0] = aold[0] - (derb2 * dera - derab * derb)/divider;
      anew[1] = aold[1] - (-derab * dera + dera2 * derb)/divider;
    }

    va = apply_funcs<square2<double>,mdiff<double>,double *,double *>(&y[0],&yhat[0],n,0);
    //va = sum_with< square2<double>, vec>(y - yhat);

    ret(j) = (n-2)*(con/va - 1);
  }
#ifdef _OPENMP
}
#endif
  }
  else{
    vec tmpX, aold(2), tmpX2, yhat, com, com2, anew(2);
    double dera, dera2, derb, derb2, derab, com3, sumcom2, divider, va;
    rowvec berow(2);
    int ij;

    for(int j = 0; j<D; ++j){
      tmpX = x.col(j);

      aold[1] = (cov(ly,tmpX)/var(tmpX))[0];
      aold[0] = my - aold[1]*sum(tmpX)/n;

      tmpX2 = tmpX % tmpX;

      yhat = exp(aold[0] + aold[1] * tmpX);
      com = y % yhat;
      com2 = yhat % yhat;
      sumcom2 = sum(com2);
      com3 = sum(com2 % tmpX);
      dera = sum(com) - sumcom2;
      dera2 = dera - sumcom2;
      derb = sum(com % tmpX) - com3;
      derb2 = sum(com % tmpX2) - 2 * sum(com2 % tmpX2);
      derab = derb - com3;

      divider = dera2 * derb2 - derab *derab;

      anew[0] = aold[0] - (derb2 * dera - derab * derb)/divider;
      anew[1] = aold[1] - (-derab * dera + dera2 * derb)/divider;

      ij =2;
      while (++ij<maxiters && std::abs(anew[0]+anew[1] - aold[0]-aold[1]) > tol ){
        aold = anew;

        yhat = exp(aold[0] + aold[1] * tmpX);
        com = y % yhat;
        com2 = yhat % yhat;
        sumcom2 = sum(com2);
        com3 = sum(com2 % tmpX);
        dera = sum(com) - sumcom2;
        dera2 = dera - sumcom2;
        derb = sum(com % tmpX) - com3;
        derb2 = sum(com % tmpX2) - 2*sum(com2 % tmpX2);
        derab = derb - com3;

        divider = dera2 * derb2 - derab *derab;

        anew[0] = aold[0] - (derb2 * dera - derab * derb)/divider;
        anew[1] = aold[1] - (-derab * dera + dera2 * derb)/divider;
      }

      //va = sum_with< square2<double>, vec>(y - yhat);
      va = apply_funcs<square2<double>,mdiff<double>,double *,double *>(&y[0],&yhat[0],n,0);
      //ret(j,0) = (n-2)*(-n * log(va) + con)/va;
      ret(j) = (n-2)*(con/va - 1);
    }
  }

  return ret;
}

mat add_term_c(const vec &y, const mat &xinc, const mat &xout, const double devi_0,
              const add_term_ini_vars& ini_vars, const double tol = 1e-07, const bool logged = false, const bool parallel = 1, const int maxiters = 100, const double ret_stat = 0){
  int nrows = xinc.n_rows;
  int selectedColumnSize = xinc.n_cols;

  int idxsz = xout.n_cols, inttype = ini_vars.inttype;
  mat out;

  if(ret_stat){
    out = mat(idxsz, 2);
  }
  else{
    out = mat(idxsz, 1);
  }

  if(parallel){
#ifdef _OPENMP
#pragma omp parallel
{
#endif
  mat tmpmat = mat(xinc);
  tmpmat.resize(nrows, selectedColumnSize+1);
  vec tmpvec;
  double out_0, out_1;
#ifdef _OPENMP
#pragma omp for
#endif
  for(int i = 0; i < idxsz; ++i){
    tmpmat.col(selectedColumnSize) = xout.col(i);
    if(inttype == 1){
      out_0 = devi_0 - glm_logistic3(tmpmat, y,&ini_vars.my[0],ini_vars.ini, tol,maxiters);

      out_1 = R::pchisq(out_0, 1, false, logged);
    }
    else if(inttype == 2){//2.0*(ylogy-sum(y%yhat));
      out_0 = devi_0 - (2*ini_vars.ylogy+glm_poisson3( tmpmat, y, ini_vars.D0, tol,maxiters));
      out_1 = R::pchisq(out_0, 1, false, logged);
    }
    else if(inttype == 3){
      // here ylogy is sum(log(y))
      out_0 = 2*(weib_reg2(y, tmpmat, ini_vars.ini, ini_vars.ylogy, tol, maxiters) - devi_0);
      out_1 = R::pchisq(out_0, 1, false, logged);
    }
    else if(inttype == 4){ // same normlog_reg
      tmpvec = prop_reg2(tmpmat, y, &ini_vars.my[0], ini_vars.ini,tol, maxiters);
      out_0 = tmpvec[2];
      out_1 = R::pchisq(out_0, 1, false, logged);
    }
    else if(inttype == 5){
      tmpvec = qpois_reg2(tmpmat, y, ini_vars.D0, ini_vars.ylogy,tol,maxiters);
      out_0 = (devi_0 - tmpvec(0))/tmpvec(1);
      out_1 = R::pf(out_0, 1, nrows-selectedColumnSize-1, false, logged);
    }
    else if(inttype == 6){ // spml
      out_0 = 2*(spml_reg2(ini_vars.u, tmpmat, tol, maxiters) - devi_0);
      out_1 = R::pchisq(out_0, 2, false, logged);
    }
    else if(inttype==7){ // multinom
      out_0 = 2 * (multinom_reg2(ini_vars.my, tmpmat, ini_vars.u,ini_vars.m0, ini_vars.b0, tol, maxiters) - devi_0);

      out_1 = R::pchisq(out_0, ini_vars.dof_mult, false, logged);
    }
    else{
      // normlog_reg
      tmpvec = normlog_reg2(y, tmpmat, ini_vars.ini, tol, maxiters);

      out_0 = (devi_0 - tmpvec(0))/tmpvec(1);
      out_1 = R::pf(out_0, 1, nrows-selectedColumnSize-1, false, logged);
    }

    if(ret_stat){
      out(i,0) = out_0;
      out(i,1) = out_1;
    }
    else{
      out(i,0) = out_1;
    }
  }
#ifdef _OPENMP
}
#endif
  }
  else{
    mat tmpmat = mat(xinc);
	tmpmat.resize(nrows, selectedColumnSize+1);
    vec tmpvec;
    double out_0, out_1;
    for(int i = 0; i < idxsz; ++i){
      tmpmat.col(selectedColumnSize) = xout.col(i);
      if(inttype == 1){
        out_0 = devi_0 - glm_logistic3(tmpmat, y,&ini_vars.my[0],ini_vars.ini, tol,maxiters);

        out_1 = R::pchisq(out_0, 1, false, logged);
      }
      else if(inttype == 2){//2.0*(ylogy-sum(y%yhat));
        out_0 = devi_0 - (2*ini_vars.ylogy+glm_poisson3( tmpmat, y, ini_vars.D0, tol,maxiters));
        out_1 = R::pchisq(out_0, 1, false, logged);
      }
      else if(inttype == 3){
        // here ylogy is sum(log(y))
        out_0 = 2*(weib_reg2(y, tmpmat, ini_vars.ini, ini_vars.ylogy, tol, maxiters) - devi_0);
        out_1 = R::pchisq(out_0, 1, false, logged);
      }
      else if(inttype == 4){ // same normlog_reg
        tmpvec = prop_reg2(tmpmat, y, &ini_vars.my[0], ini_vars.ini,tol, maxiters);
        out_0 = tmpvec[2];
        out_1 = R::pchisq(out_0, 1, false, logged);
      }
      else if(inttype == 5){
        tmpvec = qpois_reg2(tmpmat, y, ini_vars.D0, ini_vars.ylogy,tol,maxiters);
        out_0 = (devi_0 - tmpvec(0))/tmpvec(1);
        out_1 = R::pf(out_0, 1, nrows-selectedColumnSize-1, false, logged);
      }
      else if(inttype == 6){ // spml
        out_0 = 2*(spml_reg2(ini_vars.u, tmpmat, tol, maxiters) - devi_0);
        out_1 = R::pchisq(out_0, 2, false, logged);
      }
      else if(inttype==7){ // multinom
        out_0 = 2 * (multinom_reg2(ini_vars.my, tmpmat, ini_vars.u,ini_vars.m0, ini_vars.b0, tol, maxiters) - devi_0);

        out_1 = R::pchisq(out_0, ini_vars.dof_mult, false, logged);
      }
      else{
        // normlog_reg
        tmpvec = normlog_reg2(y, tmpmat, ini_vars.ini, tol, maxiters);

        out_0 = (devi_0 - tmpvec(0))/tmpvec(1);
        out_1 = R::pf(out_0, 1, nrows-selectedColumnSize-1, false, logged);
      }

      if(ret_stat){
        out(i,0) = out_0;
        out(i,1) = out_1;
      }
      else{
        out(i,0) = out_1;
      }
    }
  }
  return out;
}
