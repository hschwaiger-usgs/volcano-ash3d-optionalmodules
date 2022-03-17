#!/usr/local/bin/octave -qf
clear all;

% Load parameters from run script
nresmax = load('TC4_XY_idx.dat');
datsub  = load('TC4_XY_sub.dat');
submin=datsub(1);
submax=datsub(2);

% Over-ride here if you need to
%nresmax=3;
%submin=1;
%submax=3;

nres  = char('50000','25000','12500','06250','03125');
resdx = [0.050000 0.025000 0.012500 0.006250 0.003125];

nsub=submax-submin+1;
% Diffusion
%  subcase 1 : Explicit Diffusion in x
%        Domain = -1.0 < x < 1.0 ; -1.0 < y < 1.0 ; 0 < z < 0.025
%  subcase 2 : Explicit Diffusion in y
%        Domain = -1.0 < x < 1.0 ; -1.0 < y < 1.0 ; 0 < z < 0.025
%  subcase 3 : Explicit Diffusion in z
%        Domain = -0.25 < x < 0.25 ; -0.25 < y < 0.25 ; 0 < z < 2.0
%  subcase 4 : Crank-Nicolson in x
%        Domain = -1.0 < x < 1.0 ; -1.0 < y < 1.0 ; 0 < z < 0.025
%  subcase 5 : Crank-Nicolson in y
%        Domain = -1.0 < x < 1.0 ; -1.0 < y < 1.0 ; 0 < z < 0.025
%  subcase 6 : Crank-Nicolson in z
%        Domain = -0.25 < x < 0.25 ; -0.25 < y < 0.25 ; 0 < z < 2.0

GlobL1Errors=zeros(nsub,nresmax);
GlobL2Errors=zeros(nsub,nresmax);
GlobMCErrors=zeros(nsub,nresmax);

for subcase = submin:submax
 for inres = 1:nresmax
  dx = resdx(inres);

  ierfile1 = sprintf('DATA/TC4_XY_Sub%i_%s_err.dat',  ...
                subcase,strtrim(strtrim(nres(inres,:))));
  islfile1 = sprintf('DATA/TC4_XY_Sub%i_%s_sol.dat',  ...
                subcase,strtrim(strtrim(nres(inres,:))));

  L1L2errors = load(ierfile1);
  GlobL1Errors1(subcase,inres) = L1L2errors(1);
  GlobL2Errors1(subcase,inres) = L1L2errors(2);
  GlobMCErrors1(subcase,inres) = L1L2errors(3);
 end
end

% Write out the L1 errors for each subcase
ConvRate=zeros(submax,1);
hf = figure;
for subcase = submin:submax
 ConvRate(subcase)=log(GlobL1Errors1(subcase,1)./GlobL1Errors1(subcase,nresmax))/log(resdx(1)/resdx(nresmax));

 clf;

 x0=-1.0;
 length = 2.0;
 xm = 0.0;
 res = length/resdx(nresmax)+1;
 % Here are the material properties and initial conditions for the 1-d diffusion
 % problem.  This must be the same as in Optional_Modules/TestCases/Testcases.f90
 %     real(kind=ip) :: TC4_conc_1 = 1.0_ip  ! in kg/km3
 %     real(kind=ip) :: TC4_conc_2 = 0.0_ip  ! in kg/km3
 %     real(kind=ip) :: TC4_k_1 = 100.0_ip   ! in m2/s
 %     real(kind=ip) :: TC4_k_2 =  50.0_ip   ! in m2/s
 T1 = 1.0;
 kappa1 = 100.0*3.6e-3;
 T2 = 0.0;
 kappa2 = 50.0*3.6e-3;
 
 xx = linspace(x0,x0+length,res);
 Temperature = 0.0*xx;
 time=0.1;
 % Here's the solution from Carslaw and Jaeger 1959, p88
 Tc = (T2-T1)*(sqrt(kappa2))/(sqrt(kappa1)+sqrt(kappa2));
 for i = 1:res;
   eta1 = 0.5*(xx(i)-xm)/(sqrt(kappa1*time));
   eta2 = 0.5*(xx(i)-xm)/(sqrt(kappa2*time));
   if xx(i)<xm
     Temperature(i) = T1 + Tc*(1 + erf(eta1));
   else
     Temperature(i) = T1 + Tc*(1 + sqrt(kappa1/kappa2)*erf(eta2));
   end
 end
 dat=load(islfile1);
 x=dat(:,1);
 tsol=dat(:,2);
 csol=dat(:,3);
 err=dat(:,4);
 plot(x,tsol,'g',x,csol,'bo')
 %plot(xx-x0,Temperature,'r+')
 grid on
 legend('true','calc')

 imgfile = sprintf('PLOTS/TC4_XY_Solution_Sub%i.png',subcase);
 print (hf, imgfile, "-dpng");
end
save ("-ascii","DATA/TC4_ConvRate_XY.dat","ConvRate")

 clf;
 Ox=1.0e-1;
 Oy=1.0e-2;
 dxmin=1.0e-3;
 dxmax=1.0e-1;
 L1min=1.0e-6;
 L1max=1.0e-2;
 L2min=1.0e-6;
 L2max=1.0e-2;
 Ermin=1.0e-16;
 Ermax=1.0e-0;


 subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-2],'k--');hold on;
 subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-4],'k:')
 subplot(2,3,1),loglog(resdx(1:nresmax)',GlobL1Errors1(1:nsub,1:nresmax)');axis([dxmin dxmax L1min L1max]);hold off;
 ylabel('L1 Error')
 xlabel('dx')
 %legend('O1','O2','Exp-x','Exp-y','Exp-z','Impl-x','Impl-y','Impl-z',"location",'southeast');
 subplot(2,3,2),loglog(resdx(1:nresmax)',GlobL2Errors1(1:nsub,1:nresmax)');axis([dxmin dxmax L2min L2max]);
 ylabel('L2 Error')
 xlabel('dx')
 subplot(2,3,3),loglog(resdx(1:nresmax)',GlobMCErrors1(1:nsub,1:nresmax)');axis([dxmin dxmax Ermin Ermax]);
 ylabel('MC Error')
 xlabel('dx')

 subplot(2,3,4),loglog(resdx(1:nresmax)',GlobL1Errors1(1:nsub,1:nresmax)');axis([dxmin dxmax L1min L1max]);hold off;
 legend('Exp-x','Exp-y','Exp-z','Impl-x','Impl-y','Impl-z',"location",'southeast');

 imgfile = sprintf('PLOTS/TC4_XY.png');
 print (hf, imgfile, "-dpng");


