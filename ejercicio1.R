#Ecuaciones diferenciales parciales simultaneas
#alpha=1
#beta=2
#k12=10
#Definicion de las ecuaciones diferenciales simultaneas
#u1t=d^alpha u1/dx^alpha + k12*(u2-u1) #Ecuacion 1
#u2t=d^beta u2/dx^beta-k12*(u2-u1) #Ecuacion 1
#Definicion intervalo de las condiciones de frontera, intervalo de tiempo y los valores de  xl=0 xu=1 
#x1 < x < xu, 0 < t < tf , xl=0, xu=1
#Definicion de condiones iniciales
#Condicion inicial para u1
#u1(x, t=0)=e^(-100*(x-0.5)^2)
#Condicion inicial para u2
#u2=(x, t=0)=0
#Condicion frontera para u1
#u1(x=xl, t)=u1(x=xu, t)=0
#Condicion frontera para u2
#u2(x=xl, t)=u2(x=xu, t)=0
#Comando para borrar espacios de trabajo anteriores
rm(list = ls(all=TRUE))
#Acceso a la biblioteca para la solution numerica
library("deSolve");
setwd("~/Documentos/pde");
source("pde1a.R");
#Parametros
ncase = 1;
if(ncase==1){
  alpha=1;beta=2;k12=10;
  c_1=1;c_2=0;c_3=1;c_4=0;}
#Definicion de condiciones iniciales
f1=function(x) exp(-100*(x-0.5)^2);
f2= function(x) 0;
#Definicion de condiciones de frontera
g_0=function(t) 0;
g_L=function(t) 0;
#
#Creacion de malla espacial dx corresponde al paso discreto
xl=0; xu=1; nx=51; dx=(xu-xl)/(nx-1);
xj=seq(from=xl, to=xu, by=dx); #Dimensiones de la malla
cd=dx^(-alpha)/gamma(4-alpha);
#
#Creacion de Variable independiente para la integracion ODE paso de tiempo
t0=0;tf=0.1;nt=6;dt=(tf-t0)/(nt-1);
tout=seq(from=t0,to=tf,by=dt);
#
#Matriz A1 coeficientes a1_jk
A1=matrix(0,nrow = nx-2, ncol = nx-1 ); #Tamaño de la matriz
for(j in 1:(nx-2)){
  for (k in 0: j){
    if (k==0){
      A1[j,k+1]=(j-1)^(3-alpha)-j^(2-alpha)*(j-3+alpha);
    }else if (1<= k && k<=j-1){
      A1[j,k+1]=(j-k+1)^(3-alpha)-2*(j-k)^(3-alpha)+(j-k-1)^(3-alpha);
    }else
      A1[j, k+1]=1;
  }
}
#Matriz A2 coeficientes a1_jk
A2=matrix(0,nrow = nx-2, ncol = nx-1 ); #Tamaño de la matriz
for(j in 1:(nx-2)){
  for (k in 0: j){
    if (k==0){
      A2[j,k+1]=(j-1)^(3-beta)-j^(2-beta)*(j-3+beta);
    }else if (1<= k && k<=j-1){
      A2[j,k+1]=(j-k+1)^(3-beta)-2*(j-k)^(3-beta)+(j-k-1)^(3-beta);
    }else
      A2[j, k+1]=1;
  }
}
#Condiciones iniciales
nx=nx-2;
u0=rep(0, 2*nx);
for(j in 1:nx){
  u0[j]    = f1(xj[j+1]);
  u0[j+nx] = f2(xj[j+1]);
}
  ncall=0;
#Integracion con ODE
out = lsode(y=u0, times=tout, func=pde1a, rtol=1e-6, atol=1e-6, maxord=5);
nrow(out)
ncol(out)
#Colocacion de las matrices para u1(x, t), u2(x, t)
nx=nx+2;
u1=matrix(0, nt, nx);
u2=matrix(0, nt, nx);
#u1(x, t), u2(x, t), x ne xl xu
for(i in 1:nt){
  for(j in 2:(nx-1)){
    u1[i , j]=out[i , j];
    u2[i , j]=out[i , j+(nx-2)];
  }
}
#Reiniciar los valores de frontera
for(i in i:nt){
  u1[i, 1]=g_0(tout[i]);
  u1[i, nx]=g_L(tout[i]);
  u2[i, 1]=g_0(tout[i]);
  u2[i, nx]=g_L(tout[i]);
}
#Creacion de tabla de solucion numerica
cat(sprintf("\n\n alpha = %4.2f beta = %4.2f\n", alpha, beta));
cat(sprintf("\n   t   x   u1(x, t) u2(x, t)"));
for(i in 1:nt){
  iv=seq(from=1, to=nx, by=5);
  for(j in iv){
    cat(sprintf("\n %6.2f%6.2f%10.5f%10.5f", tout[i], xj[j], u1[i, j],u2[i, j] ))
  }
  cat(sprintf("\n"));
}
#Graficar solucion numerica de u1
matplot(xj, t(u1),type="1", lwd=2, col="black", lty=1, xlab="x", ylab="u1(x, t)", main="");
#Graficar solucion numerica de u2
matplot(xj, t(u2),type="1", lwd=2, col="black", lty=1, xlab="x", ylab="u2(x, t)", main="");
#Llamado de la funcion ODE
cat(sprintf("\n\n    ncall = %3d\n", ncall));
#