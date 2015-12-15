#' Dispersion Extraction
#' 
#' Function to return dispersion as points
#' 
#' @param q wavevector or vector of wavevectors [m^-1]
#' @param omega     angular frequency or vector of frequencies [s^-1]
#' @param eps_inf   permittivity at infinity / high frequency dielectric constant
#' @param E_f Fermi energy [eV]
#' @param tau relaxation time of electrons [s]
#' @param omega_TO vector of frequencies of transverse optical phonons [cm^-1]
#' @param omega_LO vector of frequencies of longitudinal optical phonons [cm^-1]
#' 
#' 
#' @examples 
#' q <- seq(0,55,,500)*1e6 # wavevector [m^-1]
#' omega <- seq(0,60,,500)*2*pi*1e12 # angular frequency [rad/s]
#'
#' Ef <- 0.3 # Fermi Energy
#'
#' ## Call dispersion extraction from grapheneSPP
#' ## Note: dispersion_extraction2 is a temporary function name likely to be changed in the future
#' 
#' # No phonons (actaully one down at 0, bit of a hack)
#' result<-dispersion_extraction(q,omega,E_f = Ef,omega_TO = c(0),omega_LO=c(0.1))
#' 
#' # Inclusion of SiO2 phonon frequencies [cm^-1]
#' result2<-dispersion_extraction(q,omega,E_f = Ef,
#'                                omega_TO = c(448,791.7,1128.1),omega_LO=c(498.6,811.5,1270.6))
#' 
#' ## Change units to [um^-1] and frequency [THz]
#' result$q<-result$q/1e6
#' result$omega<-result$omega/(2*pi*1e12)
#' result2$q<-result2$q/1e6
#' result2$omega<-result2$omega/(2*pi*1e12)
#' 
#' 
#' ## Plot Result
#' plot(subset(result,branch==1)$q,
#'      subset(result,branch==1)$omega,
#'      type='l',lty=2,
#'      xlab=expression(wavevector~(mu*m^-1)),
#'      ylab="Frequency (THz)",xaxs="i",yaxs="i",lwd=2)
#' 
#' lines(subset(result2,branch==1)$q,subset(result2,branch==1)$omega,type='l',col=2,lwd=2)
#' lines(subset(result2,branch==2)$q,subset(result2,branch==2)$omega,type='l',col=3,lwd=2)
#' lines(subset(result2,branch==3)$q,subset(result2,branch==3)$omega,type='l',col=4,lwd=2)
#' lines(subset(result2,branch==4)$q,subset(result2,branch==4)$omega,type='l',col=5,lwd=2)

#' @export

dispersion_extraction<-function(q,omega,eps_inf=2.4,E_f=0.3,tau=1.5e-13,omega_TO=c(448,791.7,1128.1),omega_LO=c(498.6,811.5,1270.6)){
  
  loss_matrix<-log(loss_function(q,omega,eps_inf = eps_inf,E_f = E_f,tau = tau,omega_TO = omega_TO,omega_LO = omega_LO))
  
  phonon_SO_freqs<-phonon_coupling(q,omega,d=3.5e-10,omega_TO,omega_LO,
                                   eps_graphene=3.9,
                                   eps_substrate=eps_inf)$omega_SO
  
  phonon_branch_ranges<-cut(omega,breaks = c(0,phonon_SO_freqs-1e12,Inf)) #-1e12 to correct for numerical error for branch assignment
  
  levels(phonon_branch_ranges)<-1:length(levels(phonon_branch_ranges))
  
  new_q<-c()
  new_omega<-c()
  new_branch<-c()
  
  for(bn in 1:length(levels(phonon_branch_ranges))){
    for(i in 2:length(q)){
      omega_max<-which(loss_matrix[i,]==max(loss_matrix[i,which(phonon_branch_ranges==bn)]))
      new_q<-c(new_q,q[i])
      new_omega<-c(new_omega,omega[omega_max])
      new_branch<-c(new_branch,bn)
    }
  }

result<-data.frame(q=new_q,omega=new_omega,branch=new_branch)
  return(result)
}

dispersion_extraction2<-dispersion_extraction



