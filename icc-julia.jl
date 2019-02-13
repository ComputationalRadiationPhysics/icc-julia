function parameters(Np,Nt,Ns_t)
  #parameters for calculation;
  #Np:number of the ions
  #Nt:number of the calculation circles
  #Ns_t:number of the calculation time per circle
  c=2.997*10^8;#m/s
  C=128.8;#m
  E0=931.494*10^6;#eV/u
  E=275.7*10^6;#eV/u
  gamma=1+E/E0;
  beta=sqrt(1-1/gamma^2);
  V0=beta*c;
  gamma_t=2.629;
  eta=1/gamma^2-1/gamma_t^2;
  f=beta*c/C;
  e=1.602*10^(-19);#C
  Q=5;
  A=16;
  m=A*1.66053886*10^(-27);#Kg
  Ub=10;#V
  h=1;
  T=C/V0;
  fs=(e*Q*Ub*h/(2*pi*C^2*m))^0.5;
  Tau=2.4*10^(-9);#s
  S=1;
  Lamda_opt=103.76*10^(-9);#m
  Lamda_laser=220.0055*10^(-9);#m
  Ns=Nt*Ns_t;
  delta_t=T/Ns_t;
  return c,C,E0,E,gamma,beta,V0,gamma_t,eta,f,e,Q,A,m,Ub,h,T,fs,Tau,S,Lamda_opt,Lamda_laser,Ns,delta_t
end

function initialization(Np,C,h)
  #initialize the distribution of the ion beams;
  for k=1:3
  y1=rand(Np,1)*10^(-3);
  y2=randn(Np,1)*10^3;
  for i=1:Np
    Ion[i,k]=y1[i,1]-1/2*10^(-3);
	Ion[i,3+k]=y2[i,1];
  end
  end
end

function force_3d_potential(P1_x,P1_y,P1_z,f,m)
  #calculate the bucket force as harmonic;
  #P1:position of the ion
  e=1.602*10^(-19);
  Force_bucket=Vector{Float64}(undef, 3);
  P1=[P1_x,P1_y,P1_z];
  Force_bucket=-(2*pi*f)^2*P1*m
end

function force_3d_cooling(vx,vy,vz,m,T)
  #calculate the cooling force;
  #vx,vy,vz: the velocity of the ion;
  e=2.718281828459;
  g=1/2;
  k=m/T*log(e,1/g);
  Force_cooling=Vector{Float64}(undef, 3);
  v=[vx,vy,vz];
  Force_cooling=-k*v
end
function force_3d_coulomb(P1_x,P1_y,P1_z,P2_x,P2_y,P2_z,Q1,Q2)
  #calculate the 3D coulomb force for two ions;
  #P1,P2:the position of the two ions;
  #Q1,Q2:the charge of the two ions;
  epsilon_0=8.854187817*10^(-12);
  k=1/(4*pi*epsilon_0);
  e=1.602*10^(-19);
  Force_coulomb=Vector{Float64}(undef, 3);
  P1=[P1_x,P1_y,P1_z];
  P2=[P2_x,P2_y,P2_z];
  di=P1-P2;
  disqr=(di[1]^2+di[2]^2+di[3]^2)^0.5;
  Force_coulomb=k*Q1*Q2*e^2/disqr^2*di/disqr
end

function force_3d_coulomb_all(Np,i,j,Q)
  #calculate the 3D coulomb force from all ions;
  Force_coulomb_all=Vector{Float64}(undef, 3);
  for l=1:3
    Force_coulomb_all[l]=0;
  end
  for k=1:Np
    if i!=k 
	  Force_3d_coulomb=force_3d_coulomb(Ion[(j-1)*Np+i,1],Ion[(j-1)*Np+i,2],Ion[(j-1)*Np+i,3],Ion[(j-1)*Np+k,1],Ion[(j-1)*Np+k,2],Ion[(j-1)*Np+k,3],Q,Q)
	  Force_coulomb_all=Force_coulomb_all+Force_3d_coulomb;
	else
	  Force_coulomb_all=Force_coulomb_all;
    end
  end
  return Force_coulomb_all
end

function force_3d(i,j,Q,m,T,f,Np)
  #calculate the 3D force;
  P1_x=Ion[(j-1)*Np+i,1];
  P1_y=Ion[(j-1)*Np+i,2];
  P1_z=Ion[(j-1)*Np+i,3];
  vx=Ion[(j-1)*Np+i,4];
  vy=Ion[(j-1)*Np+i,5];
  vz=Ion[(j-1)*Np+i,6];
  Force_3d_potential=force_3d_potential(P1_x,P1_y,P1_z,f,m);
  Force_3d_cooling=force_3d_cooling(vx,vy,vz,m,T);
  Force_3d_coulomb_all=force_3d_coulomb_all(Np,i,j,Q);
  Force_3d=Force_3d_potential+Force_3d_cooling+Force_3d_coulomb_all
end

function integrate_eqm(Ub,C,Q,h,Np,i,j,m,delta_t)
  pos=[Ion[(j-1)*Np+i,1],Ion[(j-1)*Np+i,2],Ion[(j-1)*Np+i,3]];
  vel=[Ion[(j-1)*Np+i,4],Ion[(j-1)*Np+i,5],Ion[(j-1)*Np+i,6]];
  calculate_acc=force_3d(i,j,Q,m,T,f,Np)/m;
  calculate_pos=pos+vel*delta_t+1/2*calculate_acc*delta_t^2;
  for l=1:3
    Ion[j*Np+i,l]=calculate_pos[l];
  end
  calculate_acc_2=force_3d(i,j+1,Q,m,T,f,Np)/m;
  calculate_vel=vel+(calculate_acc+calculate_acc_2)*delta_t*1/2;
  for l=1:3
  Ion[j*Np+i,l+3]=calculate_vel[l];
  end

end

function plotfigure(j,Np)
	fig=figure(j)
    plot(Ion[j*Np+1:j*Np+Np,1], Ion[j*Np+1:j*Np+Np,2], ".")
	xlabel("X")
	ylabel("Y")
	fig=figure(-j)
    plot(Ion[j*Np+1:j*Np+Np,1], Ion[j*Np+1:j*Np+Np,3], ".")
	xlabel("X")
	ylabel("Z")
end

Np=100;
Nt=100000;
Nout=1000;
Ns_t=10;
c,C,E0,E,gamma,beta,V0,gamma_t,eta,f,e,Q,A,m,Ub,h,T,fs,Tau,S,Lamda_opt,Lamda_laser,Ns,delta_t=parameters(Np,Nt,Ns_t);
Ion=zeros((Ns+1)*Np,6);
initialization(Np,C,h);

using PyPlot;
using ReadWriteDlm2
for j=1:Ns
for i=1:Np	
    integrate_eqm(Ub,C,Q,h,Np,i,j,m,delta_t)
end
    if rem(j,Nout)==0
      f = open( "test04_1_100circle.txt", "a" )
      for i=1:Np
        writedlm2(f,[Ion[(j-1)*Np+i,1] Ion[(j-1)*Np+i,2] Ion[(j-1)*Np+i,3] Ion[(j-1)*Np+i,4] Ion[(j-1)*Np+i,5] Ion[(j-1)*Np+i,6]]);
      end
      close( f )
	  #plotfigure(j,Np)
	end
end

#=
function plotfigure_X_Y(Np,j1,j2,delta_j)
  for j=j1:delta_j:j2
	fig=figure(j)
    plot(DD[(j-1)*Np+1:(j-1)*Np+Np,1], DD[(j-1)*Np+1:(j-1)*Np+Np,2], ".")
	xlabel("X")
	ylabel("Y")
	#xticks(-1e-3:1e-3:1e-3)
    #yticks(-1e-3:1e-3:1e-3)
  end
end
function plotfigure_X_Z(Np,j1,j2,delta_j)
  for j=j1:delta_j:j2
	fig=figure(j)
    plot(DD[(j-1)*Np+1:(j-1)*Np+Np,1], DD[(j-1)*Np+1:(j-1)*Np+Np,3], ".")
	xlabel("X")
	ylabel("Z")
	#xticks(-1e-3:1e-3:1e-3)
    #yticks(-1e-3:1e-3:1e-3)
  end
end
function closefigure(j1,j2,delta_j)
  for j=j1:delta_j:j2
  close(j);
  end
end

using PyPlot;
using ReadWriteDlm2;
DD=readdlm2("test04_1_100circle.txt")
Np=100;
j1=5;
j2=80;
delta_j=5;
plotfigure_X_Y(Np,j1,j2,delta_j)

plotfigure_X_Z(Np,j1,j2,delta_j)

closefigure(j1,j2,delta_j)

=#