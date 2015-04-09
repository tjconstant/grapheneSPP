# Internal Constants
# Reduced Planks Constant
hbar<-1.05457173e-34
# Electric Charge
e_charge<-1.60217657e-19
# Permittivity of free space
eps_0<-8.85418782e-12
# Speed of light in a vaccuum
c_0<-299792458

#' Simplifed Drude polarizability of uncoupled sheet
#' 
#' The polarizability of an uncoupled sheet described by the Drude model
#' 
#' @param q wavevector [m^-1]
#' @param omega     angular frequency [s^-1]
#' @param eps_inf   permittivity at infinity / high frequency dielectric constant
#' @param d distance between graphene and substrate [m]
#' @param E_f Fermi energy [J]
#' @param tau relaxation time of electrons [s]
#' 
#' 
#' @details
#' Note the conductivity of graphene contains two contributions for interband and intraband (Drude) transistions.
#' This model assumes h_bar*omega < 2E_f, and so interband conductivity may be neglected.
#' @return returns the total dielectirc function
#' @references
#' [1] Low, T. & Avouris, P. ACS Nano (2014), 8, 1086
#' 
#' [2] Gan, C. H. Appl. Phys. Lett. (2012) 101 111609
#' @export
pi_drude<-function(q,omega,eps_inf,E_f,tau){
  (2*(q^2)*eps_inf*E_f)/(pi^2*hbar^2*(omega^2+complex(imaginary = (omega*tau^-1))))
}

#' 2D Fourier transform of the Coulomb potential
#' 
#' 2D Fourier transform of the Coulomb potential
#' 
#' @param q wavevector [m^-1]
#' @param eps_inf   permittivity at infinity / high frequency dielectric constant
#' 
#' @return 2D Fourier transform of the Coulomb potential
#' @export
v_f<-function(q,eps_inf){
  (e_charge^2)/(2*q*eps_inf*eps_0)
}

#' Total dielectric function
#' 
#' The total dielectric function within the random phase approximation
#' 
#' @param q wavevector [m^-1]
#' @param omega     angular frequency [s^-1]
#' @param eps_inf   permittivity at infinity / high frequency dielectric constant
#' @param d distance between graphene and substrate [m]
#' @param E_f Fermi energy [J]
#' @param tau relaxation time of electrons [s]
#' 
#' @return returns the total dielectirc function
#' #' @references
#' [1] Yan, H. et al. Nat. Photonics (2013), 7, 394
#' 
#' [2] Hwang, E. H., Sensarma, R. & Das Sarma, S. Phys. Rev. B (2010), 82, 195406  
#' @export
eps_total<-function(q,omega,eps_inf,E_f,tau,omega_TO,omega_LO){
  1-v_f(q,eps_inf)*pi_drude(q,omega,eps_inf,E_f,tau)-1/(1+beta(q,omega,eps_inf = eps_inf,omega_TO = omega_TO,omega_LO = omega_LO))
}

#' Coupling Term
#' 
#' The coupling term between uncoupled conductive sheet and phonons
#' 
#' @param q wavevector [m^-1]
#' @param omega     angular frequency [s^-1]
#' @param eps_inf   permittivity at infinity / high frequency dielectric constant
#' @param d distance between graphene and substrate [m]
#' 
#' @return returns the coupling term beta
#' 
#' @references
#' [1] Wang, S Q. & Mahan G. D. Phys. Rev. B (1972) 6 4517
#' 
#' [2] Hwang, E. H., Sensarma, R. & Das Sarma, S. Phys. Rev. B (2010), 82, 195406   
#' @export
beta<-function(q,omega,d=3.5e-10,eps_inf,omega_TO,omega_LO){
  
   #omega_TO<-c(448,791.7,1128.1)
   #omega_LO<-c(498.6,811.5,1270.6)
   
   phonon_coupling(q,omega,d,omega_TO = omega_TO,omega_LO = omega_LO,eps_substrate = eps_inf)
#   
#   eps_00<-3.9
#   eps_1<-eps_00*(omega_TO[1]^2/omega_LO[1]^2)
#   eps_2<-eps_1*(omega_TO[2]^2/omega_LO[2]^2)
#   eps_3<-eps_2*(omega_TO[3]^2/omega_LO[3]^2)
#   eps_N<-2.4
#   
#   c(eps_00,eps_1,eps_2,eps_3,eps_N)
#   
#   f_1<-eps_00-eps_1
#   f_2<-eps_1-eps_2
#   f_3<-eps_2-eps_3
#   f_N<-eps_3-eps_N
#   
#   sum(c(f_1,f_2,f_3,f_N))
#   
#   omega_SO_1<-omega_LO[1]*(eps_1*(eps_00+1)/(eps_00*(eps_1+1)))^(1/2)
#   omega_SO_2<-omega_LO[2]*(eps_2*(eps_1+1)/(eps_1*(eps_2+1)))^(1/2)
#   omega_SO_3<-omega_LO[3]*(eps_3*(eps_2+1)/(eps_2*(eps_3+1)))^(1/2)
#   
#   omega_SO<-c(omega_SO_1,omega_SO_2,omega_SO_3)*c_0*100*2*pi
#   
#   alpha_1<-f_1/((eps_00+1)*(eps_1+1))
#   alpha_2<-f_2/((eps_1+1)*(eps_2+1))
#   alpha_3<-f_3/((eps_2+1)*(eps_3+1))
#   
#   return(
#     (eps_inf*exp(-2*q*d)*(
#       ((alpha_1)*(omega_SO[1])^2/(omega^2-(omega_SO[1])^2))+
#         ((alpha_2)*(omega_SO[2])^2/(omega^2-(omega_SO[2])^2))+
#         ((alpha_3)*(omega_SO[3])^2/(omega^2-(omega_SO[3])^2))
#     ))^-1)
}

#' Loss function
#'
#' Imaginary part of the inverse of the total dielectric function
#' 
#' @param q wavevector or vector of wavevectors [m^-1]
#' @param omega     angular frequency or vector of frequencies [s^-1]
#' @param eps_inf   permittivity at infinity / high frequency dielectric constant
#' @param E_f Fermi energy [eV]
#' @param tau relaxation time of electrons [s]
#' 
#' @return returns a matrix sutiable for plotting using image of the loss function as a function of both q and omega
#' @examples
#' 
#' q<-seq(0,55,,500)*1e6
#' omega<-seq(0,60,,500)*2*pi*1e12

#' image(q/1e6,omega/(2*pi*1e12),
#'      log(loss_function(q,omega,eps_inf = 2.4,E_f = 0.37,tau=1e-12)),col=rev(grey(0:100/100)),
#'      xlab=expression(list(wavevector,italic(q)~(mu*m^-1))),
#'      ylab=expression(list(frequency,italic(f)~(THz))))
#'
#'      
#' @export
loss_function <- function(q,omega,eps_inf,E_f,tau=1e-13,omega_TO=c(448,791.7,1128.1),omega_LO=c(498.6,811.5,1270.6)){
  -Im(1/outer(q,omega,eps_total,eps_inf,E_f=E_f*e_charge,tau,omega_TO,omega_LO))
}




