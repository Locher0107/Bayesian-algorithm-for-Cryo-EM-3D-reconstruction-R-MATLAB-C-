

#measure unit: nm

lambda=0.0025;#200 kv 
F=-60
Cs=1.1*10^6

x=seq(0, 8, 0.01)
ctf=as.numeric()



for(i in 1:length(x)){
  gamma_0=pi*lambda*F*x[i]^2+(1/2)*pi*Cs*lambda^3*x[i]^4
  ctf[i]=sin(gamma_0)
  
  
}
plot(x,ctf, type='l',sub='lambda=0.0025nm, C_s=1.1mm, deltaf=-60nm', xlab="frequency k: 1/nm")
title("Contrast trasfer function")



###
##Voltage= 300kv
lambda= 0.020 #A; measure unit is Angstrom
##Defocus =-2.3 mu_m
F=-23000 #A
Cs=2.7*10^7
#angstrom per pixel 1.42 A
#maximum frequency
max_f=1/2/1.42
x=seq(0, 0.35, 0.0005);
ctf=as.numeric()


for(i in 1:length(x)){
  gamma_0=pi*lambda*F*x[i]^2+(1/2)*pi*Cs*lambda^3*x[i]^4
  ctf[i]=sin(gamma_0)
  
}
plot(x, ctf, 'l')





##This is to test spyder example

###
##Voltage= 300kv
lambda= 0.020 #A; measure unit is Angstrom
##Defocus =-2.3 mu_m
F=-20000 #A
Cs=2.0*10^7
#angstrom per pixel 2.82 A
#maximum frequency
max_f=1/2/2.82
x=seq(0, 0.2, 0.0005);
ctf=as.numeric()


for(i in 1:length(x)){
  gamma_0=pi*lambda*F*x[i]^2+(1/2)*pi*Cs*lambda^3*x[i]^4
  ctf[i]=sin(gamma_0)
  
}
plot(x, ctf, 'l')






##to test coursera example

##Voltage= 300kv
lambda= 0.020 #A; measure unit is Angstrom
##Defocus =-2.3 mu_m
F=-440 #A
#F=-sqrt(4/3*Cs*lambda)
Cs=0.2*10^7
F=-sqrt(4/3*Cs*lambda)
#angstrom per pixel 2.82 A
#maximum frequency
max_f=1/2/2.82
x=seq(0, 1, 0.01);
ctf=as.numeric()


for(i in 1:length(x)){
  gamma_0=pi*lambda*F*x[i]^2+(1/2)*pi*Cs*lambda^3*x[i]^4
  ctf[i]=sin(gamma_0)
  
}
plot(x, ctf, 'l')


