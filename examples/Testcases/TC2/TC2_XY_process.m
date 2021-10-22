clear all;

% Load parameters from run script
nresmax = load('TC2_XY_idx.dat');
datsub  = load('TC2_XY_sub.dat');
submin=datsub(1);
submax=datsub(2);
datlim  = load('TC2_XY_lim.dat');
limmin=datlim(1)+1;
limmax=datlim(2)+1;

% Over-ride here if you need to
%nresmax=3;
%submin=1;
%submax=1;
%limmin=1;
%limmax=7;

nres  = char('50000','25000','12500','06250','03125');
resdx = [0.05000 0.025000 0.012500 0.006250 0.003125];

nsub=submax-submin+1;
% Vertical advection:
%  subcase 1 : Wind blows up (no fall velocity)
%  subcase 2 : Wind blows down (no fall velocity)
%  subcase 3 : No z wind (fall velocity +)
%  subcase 4 : No z wind (fall velocity -)

nlim = limmax-limmin+1;
nlim_label = char('LIM_NO','LIM_SB','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_MC');
nlimLegLabel=char('NO','SB','LW','BW','FM','MM','MC');

GlobL1Errors1=zeros(nsub,nresmax,nlim);
GlobL2Errors1=zeros(nsub,nresmax,nlim);
GlobMCErrors1=zeros(nsub,nresmax,nlim);
GlobL1Errors2=zeros(nsub,nresmax,nlim);
GlobL2Errors2=zeros(nsub,nresmax,nlim);
GlobMCErrors2=zeros(nsub,nresmax,nlim);

for subcase = submin:submax
 for ilim = limmin:limmax
  for inres = 1:nresmax
   dx = resdx(inres);
 
   step = 1;
   ierfile1 = sprintf('DATA/TC2_XY_Sub%i_%s_%s_St%i_err.dat',  ...
                 subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);
 
   L1L2errors = load(ierfile1);
   GlobL1Errors1(subcase,inres,ilim) = L1L2errors(1);
   GlobL2Errors1(subcase,inres,ilim) = L1L2errors(2);
   GlobMCErrors1(subcase,inres,ilim) = L1L2errors(3);
 
   step = 2;
   ierfile2 = sprintf('DATA/TC2_XY_Sub%i_%s_%s_St%i_err.dat',  ...
                 subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);
 
   L1L2errors = load(ierfile2);
   GlobL1Errors2(subcase,inres,ilim) = L1L2errors(1);
   GlobL2Errors2(subcase,inres,ilim) = L1L2errors(2);
   GlobMCErrors2(subcase,inres,ilim) = L1L2errors(3);
  end
 end
end

% Write out the L1 errors for each subcase
ConvRate=zeros(submax,nlim);
hf = figure;
for subcase = submin:submax
 ConvRate(subcase,:)=log(GlobL1Errors1(subcase,1,1:nlim)./GlobL1Errors1(subcase,nresmax,1:nlim))/log(resdx(1)/resdx(nresmax));

 clf;
 %Specify point of O(1) and O(2) intersection
 Ox=1.0e-1;
 Oy=1.0e-1;
 dxmin=1.0e-3;
 dxmax=1.0e-1;
 L1min=1.0e-5;
 L1max=1.0e0;
 L2min=1.0e-4;
 L2max=1.0e1;
 Ermin=1.0e-16;
 Ermax=1.0e-5;
 
 % Interior
 subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-2],'k--');hold on;
 subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-4],'k:')
 
 subplot(2,3,1),loglog(resdx(1:nresmax),GlobL1Errors1(subcase,1:nresmax,:));axis([dxmin dxmax L1min L1max]);hold off;
 ylabel('L_1 Error')
 xlabel('dx')
 %legend('O(1)','O(2)','None','Superbee','Lax-Wen','BeamWarm','Fromm','MinMod','MC','Location','northwest');
 subplot(2,3,2),loglog(resdx(1:nresmax),GlobL2Errors1(subcase,1:nresmax,:));axis([dxmin dxmax L2min L2max]);
 ylabel('L_2 Error')
 xlabel('dx')
 subplot(2,3,3),loglog(resdx(1:nresmax),GlobMCErrors1(subcase,1:nresmax,:));axis([dxmin dxmax Ermin Ermax]);
 ylabel('MC Error')
 xlabel('dx')
 
 % Boundary
 subplot(2,3,4),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-2],'k--');hold on;
 subplot(2,3,4),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-4],'k:')
 subplot(2,3,4),loglog(resdx(1:nresmax),GlobL1Errors2(subcase,1:nresmax,:));axis([dxmin dxmax L1min L1max]);hold off;
 ylabel('L_1 Error')
 xlabel('dx')
 subplot(2,3,5),loglog(resdx(1:nresmax),GlobL2Errors2(subcase,1:nresmax,:));axis([dxmin dxmax L2min L2max]);
 ylabel('L_2 Error')
 xlabel('dx')
 subplot(2,3,6),loglog(resdx(1:nresmax),GlobMCErrors2(subcase,1:nresmax,:));axis([dxmin dxmax Ermin Ermax]);
 ylabel('MC Error')
 xlabel('dx')
 
 imgfile = sprintf('PLOTS/TC2_XY_Sub%i.png',subcase);
 print (hf, imgfile, "-dpng");
end
ConvRate
save TC2_ConvRate_XY.dat ConvRate;