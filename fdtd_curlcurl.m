% This is a 1D FDTD simulation with pulse
% It displays a "movie" of the signal
% Size of the FDTD space
clear;
ke=50;
% Position of the source
ks=ke/2;
% Number of time steps
nsteps=50;
% Cell size and time stepping
c0=3.e8;
dx=0.01;
% Use Courant condition
dt=dx/(1.*c0);
% Constants
cc=c0*dt/dx;    
% Initialize vectors
Ez=zeros(1,ke);
% Hy=zeros(1,ke);
% Gaussian pulse
t0=20;
spread=8;
% Start loop
Ez_tmp = Ez;
Ez_tmp_tmp = Ez;
Ez_hist = [];

for t=1:nsteps
    
    for k=2:ke-1
        Ez(k)= 2*Ez_tmp(k)-Ez_tmp_tmp(k) + cc^2*(Ez_tmp(k+1)-2*Ez_tmp(k)+Ez_tmp(k-1));
    end
    % Source
    Ez(ks)=exp(-.5*((t-t0)/spread)^2);
    Ez_tmp_tmp = Ez_tmp;
    Ez_tmp = Ez;
    plot(Ez);axis([1 ke -2 2]);
    Ez_hist = [Ez_hist;Ez];
    xlabel('x');
    ylabel('E_z');
    pause(0.1)
end

spectrum = fft2(Ez_hist);
plot_spectrum = log(abs(spectrum));


Ft = 1/dt; 
Fx = 1/dx;

dFt = Ft/nsteps;
dFx = Fx/ke;

kx = linspace(0,Fx/2,ke/2); %Equal spaced on positive 
omega = linspace(0,Ft/2,nsteps/2); %Equal spaced in positive spectrum

kx = kx.*2*pi; %2pi
omega = omega.*2*pi; %2pi

pcolor(kx.*dx ,omega.*dx./c0 ,plot_spectrum(1:(ke/2),(1:nsteps/2)))
% pcolor(x ,y ,plot_spectrum)

% hold on
% plot(x.*dt,y.*dt,'k')

% xlabel('Frequency')
% ylabel('Phase')
shading interp
