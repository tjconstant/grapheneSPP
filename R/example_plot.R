q<-seq(0,55,,500)*1e6
omega<-seq(0,60,,500)*2*pi*1e12

image(q/1e6,omega/(2*pi*1e12),
      log(loss_function(q,omega,eps_inf = 2.4,E_f = 0.37,tau=1e-12)),col=rev(grey(0:100/100)),
      xlab=expression(list(wavevector,italic(q)~(mu*m^-1))),
      ylab=expression(list(frequency,italic(f)~(THz))))

