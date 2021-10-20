clear all;

nresmax  = 5;
nres  = char('50000', '25000','12500','06250','03125');
resdz = [0.050000 0.025000 0.012500 0.006250 0.003125];
nz    = [40 80 160 320 640];
subcase = 2;
% Vertical advection:
%  subcase 1 : vz+ vf0
%  subcase 2 : vz- vf0
%  subcase 3 : vz0 vf+
%  subcase 4 : vz0 vf-

limmin=1;
limmax=1;
nlim = limmax-limmin+1;
nlim_label = char('LIM_NO','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_SB','LIM_MC');

%for ilim = limmin:limmax
 %for inres = 1:nresmax

  ilim=1;
  inres=nresmax;

  dz = resdz(inres);
  z=linspace(0,1,nz(inres));
  step = 1;
  isolfile1 = sprintf('DATA/TC2_XY_Sub%i_%s_%s_St%i_sol.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  sol1raw = load(isolfile1);
  tsol1=sol1raw(:,1);
  csol1=sol1raw(:,2);

  step = 2;
  isolfile2 = sprintf('DATA/TC2_XY_Sub%i_%s_%s_St%i_sol.dat',  ...
                subcase,strtrim(nlim_label(ilim,:)),strtrim(nres(inres,:)),step);

  sol2raw = load(isolfile2);
  tsol2=sol2raw(:,1);
  csol2=sol2raw(:,2);

  subplot(1,2,1),plot(z,tsol1,'b-',z,csol1,'bo')
  subplot(1,2,2),plot(z,tsol2,'b-',z,csol2,'bo')

 %end
%end

