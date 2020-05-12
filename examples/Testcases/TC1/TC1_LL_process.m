clear all;

nresmax  = 5;
nres  = char('2000','1000','0500','0250','0125');
resdx = [2.0 1.0 0.5 0.25 0.125];
subcase = 8;
% Horizontal advection:
%  subcase 1 : x+y0
%  subcase 2 : x-y0
%  subcase 3 : x0y+
%  subcase 4 : x0y-
%  subcase 5 : x+y+
%  subcase 6 : x-y-
%  subcase 7 : x-y+
%  subcase 8 : x+y-

limmin=1;
limmax=7;
nlim = limmax-limmin+1;
nlim_label = char('LIM_NO','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_SB','LIM_MC');

GlobL1Errors1=zeros(nresmax,nlim);
GlobL2Errors1=zeros(nresmax,nlim);
GlobMCErrors1=zeros(nresmax,nlim);
GlobL1Errors2=zeros(nresmax,nlim);
GlobL2Errors2=zeros(nresmax,nlim);
GlobMCErrors2=zeros(nresmax,nlim);


for ilim = limmin:limmax
 for inres = 1:nresmax
  dx = resdx(inres);

  step = 1;
  ierfile1 = sprintf('DATA/TC1_LL_Sub%i_%s_%s_St%i_err.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  L1L2errors = load(ierfile1);
  GlobL1Errors1(inres,ilim) = L1L2errors(1);
  GlobL2Errors1(inres,ilim) = L1L2errors(2);
  GlobMCErrors1(inres,ilim) = L1L2errors(3);

  step = 2;
  ierfile2 = sprintf('DATA/TC1_LL_Sub%i_%s_%s_St%i_err.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  L1L2errors = load(ierfile2);
  GlobL1Errors2(inres,ilim) = L1L2errors(1);
  GlobL2Errors2(inres,ilim) = L1L2errors(2);
  GlobMCErrors2(inres,ilim) = L1L2errors(3);
 end
end

%Specify point of O(1) and O(2) intersection
Ox=1.0e1;
Oy=1.0e-4;
dxmin=1.0e-1;
dxmax=1.0e1;
L1min=1.0e-8;
L1max=1.0e-4;
L2min=1.0e-11;
L2max=1.0e-7;
Ermin=1.0e-16;
Ermax=1.0e-13;
figure;
% Interior
subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-2],'k--');hold on;
subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-4],'k:')
subplot(2,3,1),loglog(resdx(1:nresmax),GlobL1Errors1(1:nresmax,:));axis([dxmin dxmax L1min L1max]);hold off;
ylabel('L_1 Error')
xlabel('dx')
%legend('O(1)','O(2)','None','Lax-Wen','BeamWarm','Fromm','MinMod','Superbee','MC','Location','northwest');
subplot(2,3,2),loglog(resdx(1:nresmax),GlobL2Errors1(1:nresmax,:));axis([dxmin dxmax L2min L2max]);
ylabel('L_2 Error')
xlabel('dx')
subplot(2,3,3),loglog(resdx(1:nresmax),GlobMCErrors1(1:nresmax,:));axis([dxmin dxmax Ermin Ermax]);
ylabel('MC Error')
xlabel('dx')

% Boundary
subplot(2,3,4),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-2],'k--');hold on;
subplot(2,3,4),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-4],'k:')
subplot(2,3,4),loglog(resdx(1:nresmax),GlobL1Errors2(1:nresmax,:));axis([dxmin dxmax L1min L1max]);hold off;
ylabel('L_1 Error')
xlabel('dx')
subplot(2,3,5),loglog(resdx(1:nresmax),GlobL2Errors2(1:nresmax,:));axis([dxmin dxmax L2min L2max]);
ylabel('L_2 Error')
xlabel('dx')
subplot(2,3,6),loglog(resdx(1:nresmax),GlobMCErrors2(1:nresmax,:));axis([dxmin dxmax Ermin Ermax]);
ylabel('MC Error')
xlabel('dx')