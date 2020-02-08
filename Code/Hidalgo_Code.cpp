#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat rdirichlet_cpp(int num_samples,
                         arma::colvec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each row will be a draw from a Dirichlet
  arma::mat distribution = arma::zeros(distribution_size,num_samples);
  
  for (int i = 0; i < num_samples; ++i) {
    double sum_term = 0;
    // loop through the distribution and draw Gamma variables
    for (int j = 0; j < distribution_size; ++j) {
      double cur = R::rgamma(alpha_m[j],1.0);
      distribution(j,i) = cur;
      sum_term += cur;
    }
    // now normalize
    for (int j = 0; j < distribution_size; ++j) {
      distribution(j,i) = distribution(j,i)/sum_term;
    }
  }
  return(distribution);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
double log_Likelihood_double(double mu_obs, double d){
  return log(d)-(d+1)*log(mu_obs) ;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// creates a table of a vector, keeping also the absences (zeros)

// [[Rcpp::export]]
arma::vec table_factor(arma::colvec X, arma::colvec levels){
  int n = levels.size(), nX = X.size();
  arma::vec tabbozzo(n),  count(nX);
  for(int i=0; i < n; i++){
    for(int j=0; j < nX; j++){
      count[j] =  (X[j]==levels[i]);
    }
    tabbozzo[i]=accu(count);
  } 
  return(tabbozzo);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// update the distributional stick breaking weights
// [[Rcpp::export]]
arma::colvec UPD_pl(arma::colvec Ci, int L, arma::colvec alpha ,   arma::colvec lev){
  arma::colvec ul(L),  mk;
  mk   = table_factor(Ci, lev);
  arma::colvec alpha_star  = mk+alpha;
  return rdirichlet_cpp(1,alpha_star).col(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List Stratified_operations_0(arma::vec x, arma::vec col1, int val1){
  
  arma::uvec inds1      = find( col1 == val1 ) ;
  arma::vec  x10        = x.elem(inds1);
  
  int n10 = x10.size(); 
  double sl10=0.0;
  if(n10 > 0){   
    sl10 = sum(log(x10));
  }
  return List::create(
    _["sl10"]  = sl10,
    _["n10"]   = n10);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List Int_cpp(arma::colvec mu_obser, arma::vec  Ci, int L){ // performs stratified operations for every l=1,..,L
  List R;
  arma::colvec SL0(L),  nl0(L);
  for(int ls = 0; ls<L; ls++){
    R = Stratified_operations_0(mu_obser, Ci, ls+1);    
    nl0(ls) = R["n10"];
    SL0(ls) = R["sl10"];
  }
  return(List::create(
      _["SL0"]    =  SL0,
      _["nl0"]    =  nl0
  ));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double Norm_Constant_Z_l2( int Nzi_l, int N, double CSI, int q){
  int Nb = Nzi_l-1;
  int Nw = N-Nzi_l;
  double res=0; ; 
  for( int nb=0; nb<(q+1); nb++ ){
    res +=  Rf_choose(Nb, nb) * Rf_choose(Nw, q-nb) * pow(CSI, nb) * pow(1-CSI,q-nb);
  }
  return (res);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List index_row_col(arma::mat Nq, int q, int N){
  arma::umat index_row(q,N);
  arma::field<arma::uvec> index_col(N);
  for(int i=0; i<N; i++){
    index_row.col(i) = find(Nq.row(i)==1);
    arma::uvec index_col2 = find(Nq.col(i)==1);  
    index_col(i) = index_col2;
  }
  return(List::create(
      _["IR"]    =  index_row,
      _["IC"]    =  index_col
  ));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::colvec  UPD_C_Nq_MA(arma::colvec mu_obser, arma::colvec vec_all_d_l, arma::colvec pl, 
                              int L, int N,int q,arma::colvec Ci, double CSI, 
                              arma::umat index_row, arma::field<arma::uvec> index_col) {
  R_CheckUserInterrupt();
  arma::colvec possible_label = arma::linspace<arma::vec>(1, L, L);
  arma::colvec n_i_in_l(L), m_i_in_l(L), Zeta(L), corr(L), 
  log_Nq_given_z_i(L), pp(L), pp2(L), logLik_mu_i(L);
  
  
  for(int i=0; i<N; i++){  // for every label
    arma::colvec fakeC=Ci;
    logLik_mu_i.fill(0);
    n_i_in_l.fill(0);
    m_i_in_l.fill(0);
    arma:: uvec IR=index_row.col(i), 
      IC=index_col(i);
    for(int l=0; l<(L); l++){ // for every possible value taken by Ci
      logLik_mu_i(l) = log_Likelihood_double( mu_obser[i], vec_all_d_l[l]  );
      fakeC[i] = l+1;
      int Nzi_l = accu(fakeC==(l+1));
      
      //////////////////////////////////////////////////////////////////////////////////
     // int      N_index_col = index_col.n_elem;
      n_i_in_l[l] = accu( fakeC.elem(IR) == (l+1) );
      m_i_in_l[l] = accu( fakeC.elem(IC) == (l+1) );
      Zeta[l]       = Norm_Constant_Z_l2(Nzi_l,  N, CSI, q);
      corr[l]       = Norm_Constant_Z_l2(Nzi_l-1,  N, CSI, q)/Zeta[l];
      corr[l]       = (Nzi_l-1)*log(corr[l]);
      //////////////////////////////////////////////////////////////////////////////////
      //     Rcout << i << "---" << l << "---" << n_i_in_l << "\n";
      
    }
    //Rcout << n_i_in_l << "..," << m_i_in_l;
    log_Nq_given_z_i =  (n_i_in_l+m_i_in_l) * (log(CSI)-log(1-CSI)) - log(Zeta) + corr;
    pp = logLik_mu_i + log_Nq_given_z_i + log(pl);
    // closing loop on labels to compute probabilities
    if(arma::is_finite(max(pp))){
      pp2 = exp( pp - max(pp));
      Ci(i) = RcppArmadillo::sample(possible_label, 1, TRUE, pp2)[0];
    }else{
      Ci(i) = RcppArmadillo::sample(possible_label, 1, TRUE)[0];    
    }
  } // closing loop on observations
  return Ci;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REPULSIVE

// [[Rcpp::export]]
double g(double d, double t1, double t2){
  double  y = (d-t1)/t2;
  return 1 / ( 1 + exp(-y) );
  //return exp(-tau*( exp((-nu)*log(d))));
}

// [[Rcpp::export]]
double h(arma::colvec d, double t1, double t2){
  int N = d.n_rows;
  //  //Rcout << N << "***\n";
  arma::colvec dd=sort(d);
  //  //Rcout << d << "kkk\n";
  double minval=100000.;
  double st=1.;
  ////////// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
  for(int i=0;i<(N-1);i++){
    st = g( std::abs( dd(i+1) - dd(i) ),t1,t2);
    //    //Rcout << std::abs( dd(i+1) - dd(i) ) << "???\n";
    //   //Rcout << st << "--\n";
    if(st<minval){minval=st;}
  }
  
  return minval;
}

// [[Rcpp::export]]
double inverseg(double u, double t1, double t2){
  return t1 - t2 * log( 1/u - 1  ) ;
}


// [[Rcpp::export]]
double  cdf_cpp( double x, arma::mat Intervals, double a, double b){
  
  //b = 1/b; // porto da rate a scale
  int Nintervals = Intervals.n_rows;
  arma::mat    pgammaInter(Nintervals,2);
  arma::colvec deltapgamma(Nintervals+1);
  
  for(int i=0; i < Nintervals; i++){
    pgammaInter(i,0) = R::pgamma(Intervals(i,0), a, b,1,0);
    pgammaInter(i,1) = R::pgamma(Intervals(i,1), a, b,1,0); 
    if(i==0){
      deltapgamma(i) = pgammaInter(i,0);
    }else{
      deltapgamma(i) = pgammaInter(i,0)-pgammaInter(i-1,1);
    }
  }
  deltapgamma(Nintervals) = 1 - pgammaInter(Nintervals-1,1);
  arma::colvec cumprob = cumsum(deltapgamma);
  double       Cost    = accu(deltapgamma);
  // Rcout << pgammaInter<< "hhhh\n" << Cost << "''''\n" << deltapgamma << "''''\n" << cumprob;
  
  if( x <= Intervals(0,0) ){
    return( R::pgamma(x,a, b,1,0)/Cost );
  }
  
  for(int ii=0; ii<Nintervals; ii++){
    if( x> Intervals(ii,0) && x<= Intervals(ii,1) ){
      double value = cumprob(ii) / Cost;
      return(value);
    }
  }
  
  for (int iii=1; iii<Nintervals; iii++){
    if( x > Intervals(iii-1,0) && x <= Intervals(iii,1) ){
      double value = ( R::pgamma(x,a,b,1,0) - pgammaInter(iii-1,1) + cumprob[iii-1]  ) / Cost ;
      return(value);
    }
  }
  
  if( x> Intervals(Nintervals-1,1) ){
    double value = ( R::pgamma(x,a,b,1,0) - pgammaInter(Nintervals-1,1) + cumprob[Nintervals-1] ) / Cost;
    return(value);
  }
  return 0.0;
}

// [[Rcpp::export]]
double ff(double x, double u, double a_cdf, double b_cdf, arma::mat Intervals){
  //Find_the_root_of_this
  return( cdf_cpp(x, Intervals, a_cdf, b_cdf) - u);
}

// [[Rcpp::export]]
double bisection2(double u, double a_cdf, double b_cdf, arma::mat Intervals, double x0, double x1, int maxIterations=10000, double epsilon=.00001) {
  double x2=.0;
  arma::colvec init(2);
  init(0) = ff(x0,u,a_cdf,b_cdf,Intervals);
  init(1) = ff(x1,u,a_cdf,b_cdf,Intervals);
  
  if( arma::prod(sign(init))==1  ){Rcout<< "equal sign at the extremes!"; return 0.0;}
  do {
    x2 = x0*.5+x1*.5;
    init(0) = ff(x0,u,a_cdf,b_cdf,Intervals);
    init(1) = ff(x2,u,a_cdf,b_cdf,Intervals);
    if( arma::prod(sign(init))==1  ){ x0=x2;}else{x1=x2;}  
  } while (maxIterations-- > 0 && std::abs(x1 - x0) > epsilon);
  return x2;
}


// [[Rcpp::export]]
arma::mat filter_int(arma::mat limits){
  int n = limits.n_rows;
  arma::mat filtered_limits = limits;
  for(int i=1; i<n; i++){
    if ( filtered_limits(i-1,1)>filtered_limits(i,0) ){
      filtered_limits(i-1,1) = filtered_limits(i,1);
      // Rcout << "yipes!";
      limits.shed_row(i);
      filtered_limits.shed_row(i);
      i = i-1; n=n-1;
    }}
  return(filtered_limits);
}

// [[Rcpp::export]]
double Computing_Aj_sampling_dj3(int j, arma::colvec d, int L, double R, double a_post_j, 
                                 double b_post_j){
  
  arma::colvec faked = d;
  faked.shed_row(j);
  
  arma::mat limits(L-1,2);
  limits.col(0) = sort(faked)-R;
  limits.col(1) = sort(faked)+R;
  arma::mat filt_int = filter_int(limits);
  // Rcout   << limits<<"---" << a_post_j<<"---"<< 1/b_post_j<<"--\n";
  return( bisection2( R::runif(0,1), a_post_j, 1/b_post_j, filt_int, .0001,10000));
}




// [[Rcpp::export]]
arma::colvec update_d_Repulsive3(arma::colvec d, double t1, double t2, double u_rep,
                                 int L, List AIFD, double a, double b){
  
  arma::colvec faked=d;
  double R = inverseg( u_rep, t1, t2);
  arma::colvec    nl0  =  AIFD["nl0"];
  arma::colvec    SL0  =  AIFD["SL0"];
  arma::colvec    a_post = a + nl0;
  arma::colvec    b_post = b + SL0;
  for(int l=0; l<L; l++){
    // d(l) = Computing_Aj_sampling_dj2( l,  d,  L,  R, a_post(l), b_post(l));
    d(l) = Computing_Aj_sampling_dj3( l,  d,  L,  R, a_post(l), b_post(l));
  }
  return d;
}

// [[Rcpp::export]]
double update_urep(arma::colvec d,double t1, double t2){
  // << h(d,t1,t2);
  return R::runif(0,h(d,t1,t2) );
  
}  

// END REPULSIVE
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
arma::mat Tracking( arma::mat labels, arma::mat IDs, int N, int NSIM ){
  
  arma::mat trackedID(N,NSIM);
  int ind;
  for(int i =0; i<N; i++){
    for(int j=0; j<NSIM; j++){
      ind = labels(i,j);
      trackedID(i,j) =  IDs(ind-1,j);
    }}
  return trackedID;
}

// [[Rcpp::export]]
arma::colvec UPD_d( List AIFD,  double a0, double b0, int L){
  
  arma::colvec D_new(L);
  
  arma::colvec LOGSUM = as<arma::colvec>(AIFD["SL0"]);
  arma::colvec n10    = as<arma::colvec>(AIFD["nl0"]);
  arma::colvec a_star = a0 + n10;
  arma::colvec b_star = b0 + LOGSUM;
  for( int l=0; l<L; l++){
    D_new[l] = R::rgamma(a_star[l], 1/b_star[l]);
  }
  return(D_new); 
}



// [[Rcpp::export]]
double rtgamma_once(double shape, double rate, double lb, double ub) {
  // Draws from a truncated Gamma(shape, rate) distribution on (lb, ub)
  // we pass in the shape and rate parameters of the gamma
  // these are converted internally to scale parameters used by Rf_{dpq}gamma
  double lub = Rf_pgamma(lb, shape, 1.0/rate, 1,0);
  double uub = Rf_pgamma(ub, shape, 1.0/rate, 1,0);
  double u = Rf_runif(lub, uub);
  double draw = Rf_qgamma(u, shape, 1.0/rate, 1,0);
  return draw;
}

// [[Rcpp::export]]
arma::colvec UPD_d_TRUNC( List AIFD,  double a0, double b0, int L, double D){
  
  arma::colvec D_new(L);
  
  arma::colvec LOGSUM = as<arma::colvec>(AIFD["SL0"]);
  arma::colvec n10    = as<arma::colvec>(AIFD["nl0"]);
  arma::colvec a_star = a0 + n10;
  arma::colvec b_star = b0 + LOGSUM;
  for( int l=0; l<L; l++){
    D_new[l] = rtgamma_once(a_star[l], b_star[l], 0.0, D);
  }
  return(D_new); 
}

// [[Rcpp::export]]
arma::colvec UPD_d_TRUNC_MASS( List AIFD,  double a0, double b0, int L, double D, double piMass){
  
  arma::colvec D_new(L);
  arma::colvec LOGSUM = as<arma::colvec>(AIFD["SL0"]);
  arma::colvec n10    = as<arma::colvec>(AIFD["nl0"]);
  arma::colvec a_star = a0 + n10;
  arma::colvec b_star = b0 + LOGSUM;
  arma::colvec lp(2);
  for( int l=0; l<L; l++){
    double log_mx  = (::lgamma(a_star[l])) - (::lgamma(a0)) + 
                     a0*log(b0) - (a_star[l]) * log(b_star[l]) - LOGSUM[l];
    lp(0) = (log_mx)      + log(piMass)  ;
    lp(1) = log(1-piMass) + n10[l] * log(D) - (D+1)*LOGSUM[l]   ;
    lp    = exp(lp-max(lp));
    double  U = R::runif(0,1); 
    double totpeso = lp(0)/accu(lp);
    
   // Rcout << totpeso << "..." <<    "\n";
    if( U < totpeso){
      D_new[l] = rtgamma_once(a_star[l], b_star[l], 0.0, D);
    }else{
      D_new[l] = D;
    }
  }
  return(D_new); 
}


// [[Rcpp::export]]
double Norm_Constant_Z( int N, double CSI, int q, arma::colvec z, 
                        int index){
  int current_z; 
  int Nb; int Nw;
  
  int A=0 ;
  current_z = z[index];
  for( int i=0; i<N; i++ ){
    A += (z[i]==current_z);  
  }
  Nb = A-1;
  Nw = N-A;
  
  arma::colvec res(q+1); 
  for( int nb=0; nb<(q+1); nb++ ){
    //    res[nb] =  Rf_choose(Nb, nb) * Rf_choose(Nw, q-nb) * pow(CSI,nb) * pow((1-CSI),(q-nb));
    res[nb] =  Rf_choose(Nb, nb) * Rf_choose(Nw, q-nb) * pow(CSI/(1-CSI),nb);
  }
  return accu(res) * pow((1-CSI),q); // formula dopo 5.43 pag 124 tesi Facco
}


// [[Rcpp::export]]
arma::colvec nin_maker(arma::mat Nq, arma::colvec z, int N){
  arma::colvec nin(N);
  arma::rowvec a(N);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      a[j] = (z[j] == z[i]);
    }
    nin[i]       = accu(a % Nq.row(i));
  }
  return nin;
}


arma::colvec efficient_nin_maker(arma::umat index_row, arma::colvec z, int N, int q){
  arma::colvec nin(N);
  nin.fill(0);
  for(int i=0; i<N; i++){
    for(int j=0; j<q; j++){
      nin[i] = nin[i]+(z[index_row(j,i)] == z[i]);
    }
  }
  return nin;
}





// [[Rcpp::export]]
double ParetoNq_log(arma::colvec mu, arma::colvec z, arma::colvec Nin, 
                    arma::colvec d , double CSI, int q, int N){
  
  arma::colvec loglik(N);
  double current_d;
  for(int i=0;i<N;i++){
    current_d = d[  z[i]-1 ];
    loglik[i] = log(current_d) - (current_d+1.0) * log(mu[i])+
      Nin[i] * log(CSI) + (q-Nin[i]) * log(1.0-CSI) - log( Norm_Constant_Z(N, CSI, q, z,i));
  }
  return accu(loglik);
}





// [[Rcpp::export]]
List HIDALGO(arma::colvec mu_obser, 
                             int L, 
                             int NSIM, int burn_in, int thinning, bool verbose, int verbose_step, 
                             arma::mat Nq, int q, double CSI,
                             arma::colvec alpha, double a0_d, double b0_d,
                             double t1, double t2, bool REPULSIVE=0, bool TRUNC=0,bool TRUNCmass=0,  double piMass=.5, int D=0){
  
  int N = mu_obser.n_rows;         //number of observations
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  //  Initialization of Output:
  arma::mat  ALL_D(L,NSIM), ALL_Pi(L,NSIM), ALL_Ci(N,NSIM);
  ////////////////////////////////////////////////////////////////////////////////////////////
  arma::colvec  ALL_U(NSIM);
  arma::colvec  ALL_logPost(NSIM),   ALL_loglik(NSIM);
  arma::colvec  BIC(NSIM);
  arma::colvec  TAB(L);
  arma::colvec  lev  = arma::linspace<arma::vec>(1, L, L);
  ////////////////////////////////////////////////////////////////////////////////////////////
  List ICR = index_row_col(Nq,q,N);
  
  //////////////////////////////////////////////////////////////////////////////////////////
  arma::vec possible_label = arma::linspace<arma::vec>(1, L, L);
  //   Initialization of the algorithm
  arma::colvec pl = rdirichlet_cpp( 1, alpha),
    Ci = RcppArmadillo::sample(possible_label, N, TRUE, pl),
    d  = Rcpp::rgamma(L,a0_d, 1/b0_d);
  
  double u_rep = R::runif(0,1);
  
  //   ############################## GIBBS SAMPLER LOOP ################################################
  
  for(int sim=0; sim<(NSIM*thinning + burn_in) + 1; sim++){
    
    //####################################################################################################################################################################
    //# STEP 1 - sample  indicators
    //####################################################################################################################################################################for(j in 1:J){
    
    Ci = UPD_C_Nq_MA(mu_obser, d, pl,  L, N, q, Ci, CSI, ICR["IR"], ICR["IC"]);
    
    //####################################################################################################################################################################
    //# STEP 2, uptading Pi^*
    //####################################################################################################################################################################
    
    pl = UPD_pl(Ci = Ci,L = L,alpha = alpha,lev=lev);
    
    //####################################################################################################################################################################
    //# INTERVAL: updating useful quantities
    //####################################################################################################################################################################
    
    List AIFD = Int_cpp(mu_obser = mu_obser, Ci = Ci, L = L);
    
    //####################################################################################################################################################################
    //# Step 3  --  Update d
    //####################################################################################################################################################################
    
    if(TRUNC){
      d = UPD_d_TRUNC(AIFD,   a0_d,  b0_d,  L, D);
    }else  if(TRUNCmass){
      d = UPD_d_TRUNC_MASS(AIFD,   a0_d,  b0_d,  L, D, piMass);
    }else if(REPULSIVE){
      d     = update_d_Repulsive3( d,  t1,  t2,  u_rep,  L,  AIFD,  a0_d,  b0_d);//    (AIFD, a0_d, b0_d, L);
      u_rep = update_urep(d= d, t1 = t1, t2 =  t2); ////
    }else{
      d = UPD_d( AIFD,  a0_d, b0_d, L);  
    }
    
    TAB = table_factor(Ci,lev);
    
    //################################################################################################################################
    //# Storing the results
    //################################################################################################################################
    R_CheckUserInterrupt();
    
    
    if (sim > burn_in && ((sim - burn_in) % thinning == 0)) {
      int rr  = floor((sim - burn_in)/thinning);
      rr = rr-1;
      ALL_Ci.col(rr)   = Ci;
      ALL_U(rr)        = u_rep;
      ALL_Pi.col(rr)   = pl;
      ALL_D.col(rr)    = d;
      arma::colvec Nin = efficient_nin_maker(ICR["IR"], Ci, N,q);
      
      ALL_loglik(rr)   = ParetoNq_log(mu_obser, Ci, Nin, d , CSI, q, N);
      ALL_logPost(rr)  = ALL_loglik(rr) + accu(log(pl)%TAB) + sum(Rcpp::dgamma( Rcpp::NumericVector(d.begin(), d.end()), a0_d, 1./b0_d, true));
      BIC(rr)          = ALL_loglik(rr) - .5*log(N)*(N+2*L);
    }
    if (verbose) {
      if (sim % (verbose_step*thinning) == 0 && sim>0) {
        arma::colvec FF = arma::unique(Ci);
        Rcout << "Sampling iteration: "  <<  sim  << " out of " << NSIM*thinning + burn_in << "\n";
      }}
  }
  double BICm = mean(ALL_loglik) - var(ALL_loglik) * (log(N)-1);
  double AICm = mean(ALL_loglik) - var(ALL_loglik); 
  
  return List::create(
    _["AC"]  = ALL_Ci,
    _["AP"]  = ALL_Pi,
    _["AD"]  = ALL_D,
    _["AU"]  = ALL_U,
    _["AlP"] = ALL_logPost,
    _["BICm"]= BICm,
    _["AICm"]= AICm,
    _["AlLk"]= ALL_loglik,
    _["Mu"]  = mu_obser,
    _["BIC"] = BIC
  );
}
