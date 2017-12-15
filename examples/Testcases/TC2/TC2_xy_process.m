clear all;

nresmax  = 3;
nres  = char('25000','12500','06250','03125');

resdx = [0.025000 0.012500 0.006250 0.003125];
subcase = 2;
% Vertical advection:
%  subcase 1 : Wind blows up (no fall velocity)
%  subcase 2 : Wind blows down (no fall velocity)
%  subcase 3 : No z wind (fall velocity +)
%  subcase 4 : No z wind (fall velocity -)

for inres = 1:nresmax
  dx = resdx(inres);

  step = 1;
  ierfile1 = sprintf('DATA/TC2_XY_Sub%i_%s_St%i_err.dat',  ...
                subcase,strtrim(nres(inres,:)),step);

  L1L2errors = load(ierfile1);
  GlobL1Errors1(inres) = L1L2errors(1);
  GlobL2Errors1(inres) = L1L2errors(2);
  GlobMCErrors1(inres) = L1L2errors(3);

  step = 2;
  ierfile2 = sprintf('DATA/TC2_XY_Sub%i_%s_St%i_err.dat',  ...
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
subplot(1,3,3),loglog(resdx(1:nresmax),GlobMCErrors1(1:nresmax),'b-');hold on;
subplot(1,3,3),loglog(resdx(1:nresmax),GlobMCErrors2(1:nresmax),'r-');
subplot(1,3,3),loglog(resdx(1:nresmax),GlobMCErrors1(1:nresmax),'bo');
subplot(1,3,3),loglog(resdx(1:nresmax),GlobMCErrors2(1:nresmax),'ro');hold off;


ylabel('MC Error')
xlabel('dx')

