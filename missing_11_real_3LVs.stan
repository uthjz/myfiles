// model for PPMI real data, missing
// for real data
// use cholesky correlation to model random intercept and random slope
// constrain on item 1 
// item in original seq
// use 3 dimensional LVs, constrain on 3 4 11

functions {
  real expit(real x){
   real ps;
    //  xx=x[i];
      ps=exp(x)/(1+exp(x));
            return  ps;
  }


  vector convert_p(real[] psi){
    int N= size(psi); // use rows(psi) if define psi as vector
    vector [N+1] pr;
    
    pr[1]= psi[1];
   for (k in 2:N) pr[k]=psi[k]-psi[k-1];
    pr[N+1]=1-psi[N];
    return pr;
      }
  // this part is diff from full latent model
    real Sum_const_nupart (real nu_v,  real int_v, real beta_v, real cov_x_v, real u0_v) {
      //  cov_x_v as age, beta_input is len 2, u0 is len 2, nu_risk is nu len 2, int_vect is beta0
    real sum_const_nupart;
    
    sum_const_nupart=nu_v *( int_v+ cov_x_v * beta_v + u0_v); // vector operation scalar * vector
    
    return sum_const_nupart;
    
       }
       
  real Theta_Const_nupart (vector nu_vect, vector int_vect , vector beta_vect, real cov_x_v, vector u0_vect) {
      //  cov_x_v as age, beta_input is len 2, u0 is len 2, nu_risk is nu len 2, int_vect is beta0
    //int N= rows(nu_vect);
    // vector[N]  sum_const;
    real sum_const;
    sum_const=nu_vect'*( int_vect+ cov_x_v * beta_vect + u0_vect ); // vector operation scalar * vector
    
    return sum_const;
    
       }    
// following function is pointwise  value for sum(nu_k*f(s))      
   real pointV( real nu_v,  real ft_v, real tee_v, real ut_v){
     real sub_fk;
     sub_fk=nu_v* (ft_v+tee_v*ut_v);
    return sub_fk;
   }
   
   real cox_time_nupart( vector nu_vect, vector alpha_vect, vector u1_vect){
     real time_eff;
     time_eff= nu_vect'*(alpha_vect +u1_vect);
     return time_eff;
     
   }
   
   real exp_intgl(real start_t_v, real end_t_v, real time_effect_v){
     real pw_integal;
     pw_integal= (exp(time_effect_v*end_t_v)-exp(time_effect_v*start_t_v))/time_effect_v;
     return pw_integal;
     
   }
   // equal space piecewise ha
   real sum_intgl(int k_v, vector h0_vect, vector tee_sect_vect, real time_effect_v){
     // tee_sect_vect length should be k_v+1
     vector [k_v] piece_intgl;
     vector [k_v] piece_h0;
     real sum_pw;
     for(i in 1: k_v) {
              piece_intgl[i]=exp_intgl(tee_sect_vect[i], tee_sect_vect[i+1], time_effect_v);
              piece_h0[i]=h0_vect[i];
              }
              sum_pw=piece_h0'*piece_intgl;
     return  sum_pw;
     
   }
   
   vector h0_vct( int k_v, vector h0_vect){
     vector [k_v] piece_h0;
     for(i in 1: k_v) piece_h0[i]=h0_vect[i];
     
     return piece_h0;
   }
      
}
data {
  int<lower=1> num_subject;
  int<lower=1> num_obs;
//  int<lower=1> num_conti; // number of continuous outcomes
  int<lower=1> num_part1; 
  int<lower=1> num_part2; 
  int<lower=1> num_part3; 
  int<lower=1> num_ordi;  // level of each question
  int subj_long[num_obs]; // subject ID in long fmt
  
  int <lower=-1> Y_ordi_part1[num_obs, num_part1];
  int <lower=-1> Y_ordi_part2[num_obs, num_part2];
  int <lower=-1> Y_ordi_part3[num_obs, num_part3];
  
//  int conti_match[num_conti_obs] ; // if continuous obs only account part of all observation
  vector [3] a0;
 
  real<lower=0> time_obs[num_obs];  

   real age_norm[num_subject];  // normalized data
 
  real  gender_subj[num_subject];
  
  /////////////////////////////////////////////////
  // survival
  
  real <lower=0> tee [num_subject];
  int <lower=0> event[num_subject];   // possible value 0, 1, 2
  //int <lower=1, upper=num_pw> tee_id [num_subject];   // use to indicate with piece h to use

  int rr_obs [num_obs];
  // int first_last_obs [num_obs];
  int num_pw; // num of piecewise
  int <lower=1> subj_pw_ind [num_subject];//  act as indicator for pw, if t[k]<tee<t[k+1], subj_pw_ind[i]=k+1
  vector <lower=0> [num_pw+1] tee_pw [num_subject];  // 0, t1, t2...tee
  
}


parameters {
   // random intercept + random slope
  vector [6] U [num_subject];
  corr_matrix [6] Omega; // cholesky correlation matrix only intercept
  vector<lower=0> [6] Var_U; //  random scale
 // vector<lower=0> [num_conti] Var_conti;
  vector [num_obs]  ee [3];
  vector <lower=0>  [3] Var_e;
  vector [3] beta0 ; // for intercept
  vector [3] beta1 ; // for intercept
  vector [3] alpha ; // for age different to real data

//vector [num_item] a_random1; // num_theta
  vector [num_part1] a_random1; 
  vector [num_part2] a_random2;
  vector [num_part3] a_random3;


 // vector [num_item] b_random_part1  ; //use 
  vector <lower=0>  [num_part1] b_random1;  
  vector <lower=0>  [num_part2] b_random2;   
  vector <lower=0>  [num_part3] b_random3;  
  
//  vector <lower=0> [num_ordi-2] delta [num_item];
  vector<lower=0> [num_ordi-2] delta1[num_part1];
  vector<lower=0> [num_ordi-2] delta2[num_part2];
  vector<lower=0> [num_ordi-2] delta3[num_part3];

 real w;
 vector [3] eta;
 ///////////////////////////////////////////////////////////  
 // survival

  real  gam ; // survival gender num_risk 
  vector <lower=0> [num_pw]  h0 ;  // piecewise constant competing1
  vector [3] nu ;  // num_risk 
 
  
}


transformed parameters {
  cov_matrix [6] Sigma_U;
  vector<lower=0> [6] sd_U;
//  vector<lower=0>[num_conti] sd_conti;
  
//  vector [num_conti] mu_conti[num_obs]; 
  vector <lower=0> [3]  sd_e;
  ordered[num_ordi-1] a_ordi_part1[num_part1];   // this is required to use ordered_logistic
  ordered[num_ordi-1] a_ordi_part2[num_part2];
  ordered[num_ordi-1] a_ordi_part3[num_part3];
  
  
  vector   [num_part1] b_ordi_part1 ;
  vector   [num_part2] b_ordi_part2 ;
  vector   [num_part3] b_ordi_part3 ;
  

  

  vector [3] theta [num_obs];

  vector [num_part1] updrs1_ordi_hat [num_obs];
  vector [num_part2] updrs2_ordi_hat [num_obs];
  vector [num_part3] updrs3_ordi_hat [num_obs];
  
  vector [3] U0  [num_subject]; // random intercept for survival computing
  vector [3] U1  [num_subject]; // random slope for survival computing
  

  vector [num_subject] cox_const;    // constant part in cox model
  real cox_time [num_subject]  ;    // time part in cox model
  
  vector [num_subject] sum_pw_intgl;  // store sum of piece wise integral, different subject may have diff num of piece


 
/////////////////////////////////////////////////////////////
// survival computation
 real h[num_subject];
// vector  [num_pw] integral_pw1 [num_subject]; // for risk 1
// vector  [num_pw] integral_pw2 [num_subject]; // for risk 2
 
 real log_S [num_subject];
// real log_S2 [num_subject];
 real LL [num_subject];



 // nonpar & latent 
 
    for (i in 1:num_obs) {
    
     for(p in 1:3) {
     // ft[i,p] = cc[p]'*time_obs_Bt[i] ;
      U0[subj_long[i],p]= U[subj_long[i], (2*p-1)];
      U1[subj_long[i],p]= U[subj_long[i],  (2*p) ];
      theta[i,p]= beta0[p] + beta1[p]*age_norm[subj_long[i]] + alpha[p]*time_obs[i] + U0[subj_long[i],p] + U1[subj_long[i],p]*time_obs[i]+ee[p,i]; // intercept included in nonpar part
                  }
    
                       }


 //  for(i in 1: num_obs)  for (k in 1: num_conti) mu_conti[i,k]= a_conti[k]+ b_conti[k]* theta[i];
   
//////////////////////////////////////////////////////////////////////////
// ordinal 
    // set part1
 
       for (k in 1:num_part1) {
                                   a_ordi_part1[k, 1]   = a_random1[k];  
       for (lev in 2:(num_ordi-1)) a_ordi_part1[k, lev] = a_ordi_part1[k, lev-1] + delta1[k, lev-1];
                              }
    
   // set part 2
  
       for (k in 1:num_part2) {
                                   a_ordi_part2[k, 1]   = a_random2[k]; 
       for (lev in 2:(num_ordi-1)) a_ordi_part2[k, lev] = a_ordi_part2[k, lev-1] + delta2[k, lev-1];
                              }
  
  //
       for (k in 1:num_part3) {
                                   a_ordi_part3[k, 1]   = a_random3[k]; 
       for (lev in 2:(num_ordi-1)) a_ordi_part3[k, lev] = a_ordi_part3[k, lev-1] + delta3[k, lev-1];
                              }

  
// print(b_ordi_part1[1]);
// print(b_ordi_part2[1]);
    a_ordi_part1[3,1]=a0[1];
    for(lev in 2:(num_ordi-1)) a_ordi_part1[3, lev] = a_ordi_part1[3, lev-1] + delta1[3, lev-1];

    a_ordi_part2[4,1]=a0[2];
    for(lev in 2:(num_ordi-1)) a_ordi_part2[4, lev] = a_ordi_part2[4, lev-1] + delta2[4, lev-1];

    a_ordi_part3[11,1]=a0[3];
    for(lev in 2:(num_ordi-1)) a_ordi_part3[11, lev] = a_ordi_part3[11, lev-1] + delta3[11, lev-1];

    b_ordi_part1 = b_random1; // 
    b_ordi_part2 = b_random2; // 
    b_ordi_part3 = b_random3; // 

    b_ordi_part1[3]  = 1;
  
    b_ordi_part2[4]  = 1;

    b_ordi_part3[11] = 1;
  
   for (i in 1: num_obs){
   for (k in 1: num_part1) updrs1_ordi_hat[i, k]= b_ordi_part1[k]* theta[i,1];
   for (k in 1: num_part2) updrs2_ordi_hat[i, k]= b_ordi_part2[k]* theta[i,2];
   for (k in 1: num_part3) updrs3_ordi_hat[i, k]= b_ordi_part3[k]* theta[i,3];

             } 
    
///////////////////////////////////
 // survival
 
  for (i in 1:num_subject){
     
 
    
     cox_const[i]= gam*gender_subj[i] + Theta_Const_nupart(nu, beta0, beta1, gender_subj[i], U0[i]);
     cox_time[i] = cox_time_nupart(nu, alpha, U1[i]);
     
     sum_pw_intgl[i]=sum_intgl(subj_pw_ind[i], h0, tee_pw[i], cox_time[i]) ;   // integral part
     
  /////////////////////////////////////////////////////////////////////////
  // survival computation
//  pw_epart=intgl(subj_pw_ind[i], tee_pw, cox_time[i]); // not use due to not able to define length of vector
// h0_seq= h0_vct(subj_pw_ind[i], h0);
   
   
    log_S[i] =  -   exp(cox_const[i])*sum_pw_intgl[i];   //piecewise summision
      
//if (event[i]==2)      h[i] = h0[2] * exp(cox_const[i,2]+pointV(nu2[1], fb_tee[i,1], tee[i], U1[i,1]) + pointV(nu2[2], fb_tee[i,2], tee[i], U1[i,2]));// 
    if (event[i]==1)      h[i] = h0[subj_pw_ind[i]] * exp(cox_const[i]+ cox_time[i]*tee[i]); 

    if (event[i]==0)      h[i] = 1;
      
      LL[i] = log(h[i]) + log_S[i] ;  
      
                             }  
                             
//  sd_spline = sqrt(Var_spline);
  sd_e      = sqrt(Var_e);
  sd_U      = sqrt(Var_U);
  Sigma_U   = quad_form_diag(Omega, sd_U);
  
//  sd_conti= sqrt(Var_conti);

}

model {

 vector [6] zero=[0,0,0,0,0,0]';
 
 
 U ~  multi_normal(zero, Sigma_U); 
 
 
for(p in 1:3)  ee[p] ~ normal(0, sd_e[p]); // this constrain exclude other onstrains


 
  // random walk for B spline

  for(i in 1:num_obs) {
    // no missing for rr=0
    if(rr_obs[i]==0){ 
               //     for (k in 1:num_conti) Y_conti[i,k] ~ normal(mu_conti[i,k], sd_conti[k]);

                    for (k in 1: num_part1) Y_ordi_part1[i, k] ~ ordered_logistic(updrs1_ordi_hat[i, k], a_ordi_part1[k]) ;
                    for (k in 1: num_part2) Y_ordi_part2[i, k] ~ ordered_logistic(updrs2_ordi_hat[i, k], a_ordi_part2[k]) ;
                    for (k in 1: num_part3) Y_ordi_part3[i, k] ~ ordered_logistic(updrs3_ordi_hat[i, k], a_ordi_part3[k]) ;
                    
                         }

    // 
     rr_obs[i] ~ bernoulli_logit(w+ eta'* theta[i]); // evaluate all rr
    
                     }
                     
   
//  increment_log_prob(LL);
              target +=LL;
 
    beta0    ~ normal(0,20);
    beta1    ~ normal(0,20);
    
    alpha   ~ normal(0,20);
//  a_conti ~ normal(0,10);
// for(k in 1: num_conti)    b_conti [k] ~ normal(0,10);

   w ~ normal(0,10);
   eta ~ normal(0,10);
    
    for (l in 1:(num_ordi-2)) {
     //for(k in 1: num_item) delta[k,l] ~ normal(0,10) T[0,];
       for (k in 1: num_part1)   delta1[k, l] ~ normal(0, 50) T[0,] ;
       for (k in 1: num_part2)   delta2[k, l] ~ normal(0, 50) T[0,] ;
       for (k in 1: num_part3)   delta3[k, l] ~ normal(0, 50) T[0,] ;
                              }
  
 // a_random ~ normal(0,10);
  
    b_random1 ~ normal(0,20);
    b_random2 ~ normal(0, 20);
    b_random3 ~ normal(0, 20); 

    a_random1 ~ normal(0, 20);
    a_random2 ~ normal(0, 20);
    a_random3 ~ normal(0, 20);
  
  Var_U     ~ inv_gamma(0.01, 0.01);
 // Var_conti ~ inv_gamma(0.01, 0.01);
 // L_Omega~ lkj_corr_cholesky(4.0); // differe from lkj_corr()
  Omega ~ lkj_corr(2.0); //Omega=L_Omega*L_Omega'

  Var_e      ~ inv_gamma(0.01, 0.01);
 
 /////////////////////////////////////////////
 // survival
    h0 ~ gamma(0.1, 0.1); 
    nu  ~ normal(0,10);
    gam ~ normal(0,10);

}



