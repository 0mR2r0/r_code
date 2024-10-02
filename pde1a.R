  pde1a=function(t,u,parms){
#
# Function pde1a computes the derivative 
# vector of the ODEs approximating the 
# PDE
#
# Place u in two vectors, u1, u2
  u1=rep(0,nx);
  u2=rep(0,nx);
  for(i in 1:nx){
    u1[i]=u[i];
    u2[i]=u[i+nx];
  }
#
# Allocate the vectors of the ODE
# derivatives
  u1t=rep(0,nx);
  u2t=rep(0,nx);
   ut=rep(0,2*nx);
#
# Approximation of u1xx  
  u1xx=NULL;
#
# x=0
  u10=(2*dx*g_0(t)-c_2*(4*u1[1]-u1[2]))/
      (2*dx*c_1-3*c_2);
  u1xx_0=2*u10-5*u1[1]+4*u1[2]-u1[3];
  u1xx[1]=u1[2]-2*u1[1]+u10;
#
# x=1
  u1n=(2*dx*g_L(t)+c_4*(4*u1[nx]-u1[nx-1]))/
     (2*dx*c_3+3*c_4);
  u1xx[nx]=u1n-2*u1[nx]+u1[nx-1];
#
# Interior approximation of u1xx
  for(k in 2:(nx-1)){
    u1xx[k]=u1[k+1]-2*u1[k]+u1[k-1];
  }
#
# Approximation of u2xx  
  u2xx=NULL;
#
# x=0
  u20=(2*dx*g_0(t)-c_2*(4*u2[1]-u2[2]))/
      (2*dx*c_1-3*c_2);
  u2xx_0=2*u20-5*u2[1]+4*u2[2]-u2[3];
  u2xx[1]=u2[2]-2*u2[1]+u20;
#
# x=1
  u2n=(2*dx*g_L(t)+c_4*(4*u2[nx]-u2[nx-1]))/
      (2*dx*c_3+3*c_4);
  u2xx[nx]=u2n-2*u2[nx]+u2[nx-1];
#
# Interior approximation of u2xx
  for(k in 2:(nx-1)){
    u2xx[k]=u2[k+1]-2*u2[k]+u2[k-1];
  }
#
# u1 PDE
#
# Step through ODEs  
  for(j in 1:nx){
#
#   First term in series approximation of
#   fractional derivative
    u1t[j]=A1[j,1]*u1xx_0;
#
#   Subsequent terms in series approximation
#   of fractional derivative
    for(k in 1:j){
      u1t[j]=u1t[j]+A1[j,k+1]*u1xx[k];
#
#   Next k (next term in series)
    }
    u1t[j]=cd*u1t[j]+k12*(u2[j]-u1[j]);
#
# Next j (next ODE)
  }
#
# u2 PDE
#
# Step through ODEs  
  for(j in 1:nx){
#
#   First term in series approximation of
#   fractional derivative
    u2t[j]=A2[j,1]*u2xx_0;
#
#   Subsequent terms in series approximation
#   of fractional derivative
    for(k in 1:j){
      u2t[j]=u2t[j]+A2[j,k+1]*u2xx[k];
#
#   Next k (next term in series)
    }
    u2t[j]=cd*u2t[j]-k12*(u2[j]-u1[j]);
#
# Next j (next ODE)
  }
#
# Place u1t, u2t in vector ut
  for(i in 1:nx){
    ut[i]   =u1t[i];
    ut[i+nx]=u2t[i];
  }
#
# Increment calls to pde1a
  ncall <<- ncall+1;
#
# Return derivative vector of ODEs
  return(list(c(ut)));
  }