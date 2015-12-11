#' Dispersion Extraction
#' 
#' Function to return sdispersion as points
#' 
#' @param q wavevector [m^-1]
#' @param omega     angular frequency [s^-1]
#' @param loss_matrix the log of the loss matrix
#' @note Currently only 3 phonons supported
#' 
#' @examples 
#' omega<-seq(0,60,,500)*2*pi*1e12
#' q<-seq(0,2*pi/550e-9,,500)

#' plot_it<-function(Ef){
#'   loss_matrix<-log(loss_function(q,omega,eps_inf = 2.4,E_f = Ef,tau=1.5e-13))
#'   result<-dispersion_extraction(q,omega,loss_matrix)
#'   ggplot(result,aes(q,omega,color=factor(branch)))+geom_line(size=2)+theme_bw()+ylim(range(omega))
#' }
#' library(manipulate)
#' manipulate(plot_it(E_f),E_f=slider(min = 0.000001,max = 1))

#' @export

dispersion_extraction2<-function(q,omega,eps_inf=2.4,E_f=0.3,tau=1.5e-13,omega_TO=c(448,791.7,1128.1),omega_LO=c(498.6,811.5,1270.6)){
  
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



