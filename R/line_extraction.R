#' Dispersion Extraction
#' 
#' Function to return sdispersion as points
#' 
#' @param q wavevector [m^-1]
#' @param omega     angular frequency [s^-1]
#' @param loss_matrix the log of the loss matrix
#' @note Currently only 3 phonons supported
#' @export

dispersion_extraction<-function(q,omega,loss_matrix){
  
  maxima<-data.frame()
  for(i in 2:length(omega)){
  
   q_max<-q[which(loss_matrix[,i]==max(loss_matrix[,i],na.rm = T))]
   
   if(i>3){

   } 
   
   if(max(loss_matrix[,i],na.rm = T)>0){
     maxima<-rbind(maxima,
                   data.frame(
                     q=q_max,omega=omega[i]))
   }

  }
  
  return(maxima)
}





