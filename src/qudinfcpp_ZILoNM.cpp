#include <Rcpp.h>
#include <iostream>
#include <math.h>
using namespace Rcpp;

/* old: without beta5
 double g1(double m, double x, double y, double mu, double sig,
 double ypart, double beta1, double delta, double eta){
 return pow(eta,4)*delta*delta/(2*beta1*beta1) + eta*eta*(ypart-y)/beta1 - pow(beta1*m + (eta*eta*delta*delta)/beta1 + ypart-y,2)/(2*delta*delta) ;}

 double g2(double m, double mu, double sig){
 return -pow(log(m)+sig*sig - mu,2)/(2*sig*sig) + sig*sig/2.0 - mu    ;}

 double g(double m, double x, double y, double mu, double sig,
 double ypart, double beta1, double delta, double eta){
 return g1(m,x,y,mu,sig,ypart,beta1,delta,eta)+g2(m,mu,sig) ;}

 double f(double m, double x, double y, double mu, double sig,
 double ypart, double beta1, double delta, double eta, double g_star){
 return exp(g(m,x,y,mu,sig,ypart,beta1,delta,eta)- g_star) ;}
 */

double g1_original(double m, double x, double y, double mu, double sig,
                   double ypart, double beta1, double beta5, double delta, double eta){
  return - pow(y-beta1*m-beta5*x*m-ypart,2)/(2*delta*delta) - eta*eta*m;}

double g2_original(double m, double mu, double sig){
  return -log(m) - pow(log(m) - mu,2)/(2*sig*sig);}

double g_original(double m, double x, double y, double mu, double sig,
                  double ypart, double beta1, double beta5, double delta, double eta){
  return  g1_original(m,x,y,mu,sig,ypart,beta1,beta5,delta,eta)+g2_original(m,mu,sig);}

double f_original(double m, double x, double y, double mu, double sig,
                  double ypart, double beta1, double beta5, double delta, double eta, double g_star){
  return exp(g_original(m,x,y,mu,sig,ypart,beta1,beta5,delta,eta)- g_star) ;}

double fquadinf(double yy, double xa, double xb,
                double x, double y, double mu, double sig,
                double ypart, double beta1, double beta5, double delta, double eta, double g_star) {
  double u = -1;
  double v = 1 ;
  double duv = (xb-xa)/(v-u);

  return duv * f_original(xa + duv*(yy-u), x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star);
}

double quadinfcpp(double xa, double xb, List nodes, List weights,
                  double x, double y, double mu, double sig,
                  double ypart, double beta1, double beta5, double delta, double eta, double g_star) {
  double Q;
  double h=0.5;
  NumericVector n=nodes[0];
  NumericVector w=weights[0];
  double s= w[6] * fquadinf(n[6],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star);
  for (int j=0;j<6;j++) {
    s=s + w[j] * (fquadinf(n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) +
      fquadinf(-n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) );
  }
  Q=s*h;

  double delta1;
  double newQ;
  double tol=1e-12;
  for (int k=1;k<7;k++) {
    n=nodes[k];
    w=weights[k];
    s=0;
    for (int j=0;j<w.size();j++) {
      s=s+w[j]*(fquadinf(n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) +
        fquadinf(-n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) );
    }
    h = h/2;
    newQ=s*h + Q/2.0;
    delta1 = std::abs(newQ-Q);
    Q=newQ;
    if (delta1<tol) {
      break;
    }
  }

  return Q;}


double loghicpp(double x, double y, double mu, double sig,
                double ypart, double beta1, double beta5, double delta, double eta,
                List nodes, List weights, double C){
  double m1 = std::max( -((eta*eta*delta*delta)/(beta1+beta5*x) + ypart-y)/(beta1+beta5*x)  ,0.0);
  double m2 = exp(-sig*sig + mu);
  double g2_star = g2_original(m2,mu,sig);
  double start=0;
  double lower=0, upper=0;
  double m3_1=0, m3_2=0, m4_1=0, m4_2=0, g1_star=0,g_star=0,val=0,by=0;
  NumericVector gvalue(1000);

  // situation 1
  if((m1 <= C) and (m2 <= C)){
    if(m1 == 0){start=1e-200;
    }else{start=m1;}
    val=std::min(start,m2);
    by = fabs(m2 - start)/999.0;
    gvalue[0] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;
    for (double l=1;l<1000;l++){
      val = val+by;
      gvalue[l] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;}
    g_star = max(gvalue);

    //int count = 0;
    double g2_m1=g2_original(m1,mu,sig);
    double g1_m2=g1_original(m2,x,y,mu,sig,ypart,beta1,beta5,delta,eta);

    if(m1==0){
      g1_star = g1_original(0,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
      m3_1 = 0;
      if((y-ypart)/(beta1+beta5*x) <= 0){
        m4_1 = (30+g2_star-g_star-pow(y-ypart,2)/(2*delta*delta) )/(eta*eta);
      }else{
        m4_1 = (30+g2_star-g_star)/(eta*eta);
      }

    }else{
      g1_star = g1_original(m1,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
      NumericVector g130 = NumericVector::create(
        (-sqrt(2*delta*delta)*sqrt(pow(eta,4)*delta*delta/(2*pow(beta1+beta5*x,2))+eta*eta*(ypart-y)/(beta1+beta5*x)-g_star+30+g2_star) -
          (eta*eta*delta*delta/(beta1+beta5*x) +ypart-y) )/(beta1+beta5*x),
          (sqrt(2*delta*delta)*sqrt(pow(eta,4)*delta*delta/(2*pow(beta1+beta5*x,2))+eta*eta*(ypart-y)/(beta1+beta5*x)-g_star+30+g2_star) -
            (eta*eta*delta*delta/(beta1+beta5*x) +ypart-y) )/(beta1+beta5*x) );

      m3_1 = min(g130);
      m4_1 = max(g130);
    }
    //smaller roots for g1+g2*-g*=-30, g2+g1*-g*=-30 for lower bound
    //larger roots for g1+g2*-g*=-30, g2+g1*-g*=-30 for upper bound
    NumericVector expo30 = NumericVector::create(-sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu ,
                                                 sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu );
    m3_2 = exp(min(expo30));
    m4_2 = exp(max(expo30));

    // m1 < m2
    if (m1 < m2){
      // lower bound
      if (m1==0){
        lower = m3_2 ;
      }else if( g2_m1+g1_star-g_star < -30){ //check m1=m3
        lower = m1;
      }else{
        lower = std::max(m3_1,m3_2);
      }
      // upper bound
      if (g1_m2+g2_star-g_star < -30){//check m2=m4
        upper = m2;
      }else{
        upper = std::min(m4_1,m4_2);
      }
    }else{// m1 > m2
      // count = count +1;
      // lower bonund
      if(g1_m2+g2_star-g_star < -30){ //check m2=m4
        lower = m2;
      }else{
        lower = std::max(m3_1,m3_2);
      }
      // upper bound
      if (g2_m1+g1_star-g_star < -30){//check m1=m3
        upper= m1;
      }else{
        upper = std::min(m4_1,m4_2);
      }

    }
    upper = std::min(upper,C);
    // situation 2
  }else if((m1>C) != (m2>C)){
    start = std::min(m1,m2);
    if(start == 0){start=1e-200;}
    val=start;
    by = fabs(C - start)/999.0;
    gvalue[0] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;
    for (double l=1;l<1000;l++){
      val = val+by;
      gvalue[l] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;}
    g_star = max(gvalue);

    upper = C;
    if(m1<C){g1_star=g1_original(m1,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    }else{g1_star=g1_original(C,x,y,mu,sig,ypart,beta1,beta5,delta,eta);}
    lower = exp(-sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu);
    // situation 3
  }else{
    g_star = g_original(C,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    g1_star = g1_original(C,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    upper = C;
    lower = exp(-sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu);
  }
  //Rcout << "m1=" << m1;
  //Rcout << " m2=" << m2;
  //Rcout << "lower=" << lower;
  //Rcout << " upper=" << upper <<"\n";

  double integ =  quadinfcpp(lower,upper,nodes,weights,
                             x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star);
  //Rcout << "lower=" << lower;
  //Rcout << " upper=" << upper;
  //Rcout << " gs=" << g_star;
  //Rcout << " integ=" << integ <<"\n";
  double output = log(integ) + g_star;

  return output;
}

/* old: mixture with two components
 List loghicpp_all2(NumericVector X_group2, NumericVector Y_group2,
 NumericVector mu1, NumericVector sig1,NumericVector mu2, NumericVector sig2,
 NumericVector Ypart, double beta1, double delta, double eta,
 List nodes, List weights, double C){
 R_xlen_t len = X_group2.size() ;
 NumericVector loghi1(len);
 NumericVector loghi2(len);

 for(int i=0;i<len;i++){
 loghi1[i] = loghicpp(X_group2[i], Y_group2[i], mu1[i], sig1[i], Ypart[i], beta1, delta, eta, nodes, weights,C);
 loghi2[i] = loghicpp(X_group2[i], Y_group2[i], mu2[i], sig2[i], Ypart[i], beta1, delta, eta, nodes, weights,C);
 }
 List L = List::create(Named("loghi1")=loghi1, _["loghi2"] = loghi2);

 return L;
 }
 */

// [[Rcpp::export]]
NumericMatrix loghicpp_all(NumericVector X_group2, NumericVector Y_group2,
                           NumericMatrix mu, double sig,
                           NumericVector Ypart, double beta1, double beta5, double delta, double eta,
                           List nodes, List weights, double C){
  NumericMatrix loghik = clone(mu);

  for(int k=0;k < mu.ncol();k++){
    for(int i=0;i < mu.nrow();i++){
      loghik(i,k) = loghicpp(X_group2[i], Y_group2[i], mu(i,k), sig, Ypart[i], beta1, beta5, delta, eta, nodes, weights,C);
    }
  }

  return loghik;
}
