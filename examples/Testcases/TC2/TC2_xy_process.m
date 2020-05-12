clear all;

nresmax  = 4;
nres  = char('25000','12500','06250','03125');

resdx = [0.025000 0.012500 0.006250 0.003125];
subcase = 3;

nlim = 7;
nlim_label = char('LIM_NO','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_SB','LIM_MC');

% Vertical advection:
%  subcase 1 : Wind blows up (no fall velocity)
%  subcase 2 : Wind blows down (no fall velocity)
%  subcase 3 : No z wind (fall velocity +)
%  subcase 4 : No z wind (fall velocity -)

GlobL1Errors1=zeros(nresmax,nlim);
GlobL2Errors1=zeros(nresmax,nlim);
GlobMCErrors1=zeros(nresmax,nlim);
GlobL1Errors2=zeros(nresmax,nlim);
GlobL2Errors2=zeros(nresmax,nlim);
GlobMCErrors2=zeros(nresmax,nlim);

for ilim = 1:nlim
 for inres = 1:nresmax
  dx = resdx(inres);

  step = 1;
  ierfile1 = sprintf('DATA/TC2_XY_Sub%i_%s_%s_St%i_err.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  L1L2errors = load(ierfile1);
  GlobL1Errors1(inres,ilim) = L1L2errors(1);
  GlobL2Errors1(inres,ilim) = L1L2errors(2);
  GlobMCErrors1(inres,ilim) = L1L2errors(3);

  step = 2;
  ierfile2 = sprintf('DATA/TC2_XY_Sub%i_%s_%s_St%i_err.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  L1L2errors = load(ierfile2);
  GlobL1Errors2(inres,ilim) = L1L2errors(1);
  GlobL2Errors2(inres,ilim) = L1L2errors(2);
  GlobMCErrors2(inres,ilim) = L1L2errors(3);
 end
end

dxmin=1.0e-3;
dxmax=1.0e-1;
L1min=1.0e-5;
L1max=1.0e-1;
L2min=1.0e-4;
L2max=1.0e0;
Ermin=1.0e-16;
Ermax=1.0e-13;
figure;
subplot(2,3,1),loglog(resdx(1:nresmax),GlobL1Errors1(1:nresmax,:));%axis([dxmin dxmax L1min L1max]);
ylabel('L1 Error')
xlabel('dx')
%legend('None','Lax-Wen','BeamWarm','Fromm','MinMod','Superbee','MC');
subplot(2,3,2),loglog(resdx(1:nresmax),GlobL2Errors1(1:nresmax,:));%axis([dxmin dxmax L2min L2max]);
ylabel('L2 Error')
xlabel('dx')
subplot(2,3,3),loglog(resdx(1:nresmax),GlobMCErrors1(1:nresmax,:));%axis([dxmin dxmax Ermin Ermax]);
ylabel('MC Error')
xlabel('dx')

subplot(2,3,4),loglog(resdx(1:nresmax),GlobL1Errors2(1:nresmax,:));%axis([dxmin dxmax L1min L1max]);
ylabel('L1 Error')
xlabel('dx')
subplot(2,3,5),loglog(resdx(1:nresmax),GlobL2Errors2(1:nresmax,:));%axis([dxmin dxmax L2min L2max]);
ylabel('L2 Error')
xlabel('dx')
subplot(2,3,6),loglog(resdx(1:nresmax),GlobMCErrors2(1:nresmax,:));%axis([dxmin dxmax Ermin Ermax]);
ylabel('MC Error')
xlabel('dx')

