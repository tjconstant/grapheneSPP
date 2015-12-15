#' Phonon Coupling Position and Strength Calculation
#' 
#' Function to return surface phonon frequencies and coupling coupling coefficents for dispersion
#' 
#' @param q wavevector [m^-1]
#' @param omega     angular frequency [s^-1]
#' @param d distance between graphene and substrate [m]
#' @param omega_TO vector of frequencies of transverse optical phonons [cm^-1]
#' @param omega_LO vector of frequencies of longitudinal optical phonons [cm^-1]
#' @param eps_graphene  static permittivity 
#' @param eps_substrate  permittivity at infinity / high frequency dielectric constant (also refered to elsewhere as eps_inf)
#' 
#' @export
phonon_coupling<-function(q,omega,d=3.5e-10,omega_TO,
                          omega_LO,
                          eps_graphene=3.9,
                          eps_substrate=2.4){
  
  eps_phonons<-c(eps_graphene)
  for(n in 1:length(omega_TO)){
    eps_phonons<-c(eps_phonons,eps_phonons[n]*(omega_TO[n]^2/omega_LO[n]^2))
  }
  eps_phonons<-c(eps_phonons,eps_substrate)
  
  f_phonons<-c()
  for(i in 1:(length(eps_phonons)-1)){
    f_phonons<-c(f_phonons,
                 eps_phonons[i]-eps_phonons[i+1])
  }
  
  omega_SO<-c()
  alpha_phonons<-c()
  
  for(i in 1:length(omega_TO)){
    omega_SO<-c(omega_SO,
    omega_LO[i]*(eps_phonons[i+1]*(eps_phonons[i]+1)/(eps_phonons[i]*(eps_phonons[i+1]+1)))^(1/2)
    )
    
    alpha_phonons<-c(alpha_phonons,
             f_phonons[i]/((eps_phonons[i]+1)*(eps_phonons[i+1]+1))
             )
  }
  
  
  omega_SO<-omega_SO*c_0*100*2*pi
  
  #sum_vector<-c()
  #for(i in 1:length(omega_TO)){
  #               sum_vector<-c(sum_vector,
  #                             alpha_phonons[i]*omega_SO[i]^2/(omega^2-omega[i]^2))
  #}
  
  #sum_matrix<-matrix(sum_vector,nrow=length(omega_TO))
  
  #return((eps_substrate*exp(-2*q*d)*rowSums(sum_matrix))^-1)
  
  phonon_sum<-rep(0,length(omega))
  for(i in 1:length(omega_SO)){
    phonon_sum<-phonon_sum+((alpha_phonons[i])*(omega_SO[i])^2/(omega^2-(omega_SO[i])^2))
  }
  
  # This used to be the phonon sum, left here in case if ever needs to go back in
  # ((alpha_phonons[1])*(omega_SO[1])^2/(omega^2-(omega_SO[1])^2))+
  #  ((alpha_phonons[2])*(omega_SO[2])^2/(omega^2-(omega_SO[2])^2))+
  #  ((alpha_phonons[3])*(omega_SO[3])^2/(omega^2-(omega_SO[3])^2))
  
  return(
    list(
      beta=
        (eps_substrate*exp(-2*q*d)*(
          phonon_sum
          ))^-1,
      omega_SO=omega_SO)
  )
    
}

