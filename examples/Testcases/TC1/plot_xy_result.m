clear all;

nresmax  = 3;
nres  = char('50000','25000','12500','06250','03125');
resdx = [0.05000 0.025000 0.012500 0.006250 0.003125];
nx    = [40 80 160 320 640];
ny    = [40 80 160 320 640];
subcase = 5;
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
limmax=1;
nlim = limmax-limmin+1;
nlim_label = char('LIM_NO','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_SB','LIM_MC');

%for ilim = limmin:limmax
 %for inres = 1:nresmax

  ilim=7;
  inres=nresmax;

  dx = resdx(inres);
  x=linspace(-1,1,nx(inres));
  y=linspace(-1,1,ny(inres));
  step = 1;
  isolfile1 = sprintf('DATA/TC1_XY_Sub%i_%s_%s_St%i_sol.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  sol1raw = load(isolfile1);
  dat1=reshape(sol1raw(:,2),nx(inres),ny(inres));

  step = 2;
  isolfile2 = sprintf('DATA/TC1_XY_Sub%i_%s_%s_St%i_sol.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  sol2raw = load(isolfile2);
  dat2=reshape(sol2raw(:,2),nx(inres),ny(inres));
  subplot(1,2,1),contour(x,y,dat1);axis square
  subplot(1,2,2),contour(x,y,dat2);axis square

 %end
%end

