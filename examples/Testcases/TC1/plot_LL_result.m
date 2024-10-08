clear all;

nresmax  = 3;
nres  = char('2000','1000','0500','0250','0125');
resdx = [2.0 1.0 0.5 0.25 0.125];
nx    = [31 72 144 288 572];
ny    = [31 72 144 288 572];

subcase = 1;
% Horizontal advection:
%  subcase 1 : x+y0
%  subcase 2 : x-y0
%  subcase 3 : x0y+
%  subcase 4 : x0y-
%  subcase 5 : x+y+
%  subcase 6 : x-y-
%  subcase 7 : x-y+
%  subcase 8 : x+y-

nlim_label = char('LIM_NO','LIM_SB','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_MC');

ilim=1;
inres=nresmax;

dx = resdx(inres);
x=linspace(324,416,nx(inres));
y=linspace(-36,36,ny(inres));
step = 1;
isolfile1 = sprintf('DATA/TC1_LL_Sub%i_%s_%s_St%i_sol.dat',  ...
              subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);
sol1raw = load(isolfile1);
dat1=reshape(sol1raw(:,2),nx(inres),ny(inres));
tru1=reshape(sol1raw(:,1),nx(inres),ny(inres));

step = 2;
isolfile2 = sprintf('DATA/TC1_LL_Sub%i_%s_%s_St%i_sol.dat',  ...
              subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);
sol2raw = load(isolfile2);
dat2=reshape(sol2raw(:,2),nx(inres),ny(inres));
tru2=reshape(sol2raw(:,1),nx(inres),ny(inres));

subplot(1,2,1),contour(x,y,tru1);axis square
subplot(1,2,2),contour(x,y,tru2);axis square
