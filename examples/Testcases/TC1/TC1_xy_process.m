clear all;

nresmax  = 3;
nres  = char('50000','25000','12500','06250','03125');

resdx = [0.05000 0.025000 0.012500 0.006250 0.003125];
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

for inres = 1:nresmax
  dx = resdx(inres);

  step = 1;
  ierfile1 = sprintf('DATA/TC1_XY_Sub%i_%s_St%i_err.dat',  ...
                subcase,strtrim(nres(inres,:)),step);

  L1L2errors = load(ierfile1);
  GlobL1Errors1(inres) = L1L2errors(1);
  GlobL2Errors1(inres) = L1L2errors(2);
  GlobMCErrors1(inres) = L1L2errors(3);

  step = 2;
  ierfile2 = sprintf('DATA/TC1_XY_Sub%i_%s_St%i_err.dat',  ...
                subcase,strtrim(nres(inres,:)),step);

  L1L2errors = load(ierfile2);
  GlobL1Errors2(inres) = L1L2errors(1);
  GlobL2Errors2(inres) = L1L2errors(2);
  GlobMCErrors2(inres) = L1L2errors(3);

end

figure;
subplot(1,3,1),loglog(resdx(1:nresmax),GlobL1Errors1(1:nresmax),'b-',resdx(1:nresmax),GlobL1Errors2(1:nresmax),'r-')
ylabel('L1 Error')
xlabel('dx')
legend('Stp1-Inter','Stp2-Bound')
subplot(1,3,2),loglog(resdx(1:nresmax),GlobL2Errors1(1:nresmax),'b-',resdx(1:nresmax),GlobL2Errors2(1:nresmax),'r-')
ylabel('L2 Error')
xlabel('dx')
subplot(1,3,3),loglog(resdx(1:nresmax),GlobMCErrors1(1:nresmax),'b-',resdx(1:nresmax),GlobMCErrors1(1:nresmax),'bo');hold on;
subplot(1,3,3),loglog(resdx(1:nresmax),GlobMCErrors2(1:nresmax),'r-',resdx(1:nresmax),GlobMCErrors2(1:nresmax),'ro');hold off;
ylabel('MC Error')
xlabel('dx')

