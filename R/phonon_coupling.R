phonon_coupling<-function(q,omega,d=3.5e-10,omega_TO=c(448,791.7,1128.1),
                          omega_LO=c(498.6,811.5,1317),
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
  
  for(i in 1:length(omega_TO)){
                 alpha_phonons[i]*omega_SO[i]^2/(omega^2-omega[i]^2)
  }

  return(
    (eps_substrate*exp(-2*q*d)*(
      ((alpha_phonons[1])*(omega_SO[1])^2/(omega^2-(omega_SO[1])^2))+
        ((alpha_phonons[2])*(omega_SO[2])^2/(omega^2-(omega_SO[2])^2))+
        ((alpha_phonons[3])*(omega_SO[3])^2/(omega^2-(omega_SO[3])^2))
    ))^-1)
  
}

