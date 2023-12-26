#!/usr/bin/octave -qf
clear all;

IsOct    = 0;
TestCase = 6;

% Load parameters from run script
nresmax = load('TC6_XY_idx.dat');
datlim  = load('TC6_XY_lim.dat');
limmin=datlim(1)+1;
limmax=datlim(2)+1;

% Set this to 1 if you want plots of the solution and errors.
% Set to 0 if you only want convergence plots and time plots.
cploterr1 = 1; % Plot L1 error
cploterr2 = 0; % Plot L2 error
cplotmass = 1; % Plot Mass Consv error

resdx = [20.00 10.00 5.00 2.50 1.25];
resnx = [10 20 40 80 160];
resny = [10 20 40 80 160];

nlim = limmax-limmin+1;
nlim_label = char('LIM_NO','LIM_SB','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_MC');
nlimLegLabel=char('NO','SB','LW','BW','FM','MM','MC');
nres  = char('200000','100000','050000','025000','012500');

xmin = -1.0;
xmax =  1.0;
ymin = -1.0;
ymax =  1.0;

GlobL1Errors  = zeros(nresmax,nlim);
GlobL2Errors  = zeros(nresmax,nlim);
GlobMCErrors  = zeros(nresmax,nlim);
GlobExecTimes = zeros(nresmax,nlim);

cinter = linspace(0.05,0.95,10);

for il = 1:nlim
  for inres = 1:nresmax

    ierfile = sprintf('DATA/TC%i_XY_%s_%s_err.dat',  ...
        TestCase,nlim_label(il,:),nres(inres,:));
    %itimefile = sprintf('DATA/TC%i_XY_%s_etime.dat',  ...
    %    TestCase,nlim_label(il,:));

    L1L2errors = load(ierfile);
    GlobL1Errors(inres,il) = L1L2errors(1);
    GlobL2Errors(inres,il) = L1L2errors(2);
    GlobMCErrors(inres,il) = L1L2errors(3);

  end
end

% Write out the L1 errors for each limiter
ConvRate=zeros(nlim,1);
ConvRate(1:nlim)=log(GlobL1Errors(1,1:nlim)./GlobL1Errors(nresmax,1:nlim))/log(resdx(1)/resdx(nresmax));

save ("-ascii","DATA/TC6_ConvRate_XY.dat","ConvRate")

plotflags = [cploterr1, cploterr2, cplotmass];
for iplot = 1:3
  if plotflags(iplot)==1
    if iplot==1
      plotvar = GlobL1Errors;
      title_str = 'L_1 Error';
      out_file = 'PLOTS/TC5_XY_Solution_Error_L1.png';
      paxis = [1.0e0 1.0e2 1.0e-5 1.0e1];
    elseif iplot==2
      plotvar = GlobL2Errors;
      title_str = 'L_2 Error';
      out_file = 'PLOTS/TC5_XY_Solution_Error_L2.png';
      paxis = [1.0e-3 1.0e-1 1.0e-2 1.0e0];
    elseif iplot==3
      plotvar = GlobMCErrors;
      title_str = 'Mass Cons. Ratio';
      out_file = 'PLOTS/TC5_XY_Solution_Error_MassCons.png';
      paxis = [1.0e0 1.0e2 1.0e-16 1.0e1];
    end
    figure;
    subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,1,1),'-+r')
    hold on;
    if nlim > 1
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,2,1),'-og')
    end
    if nlim > 2
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,3,1),'-^b')
    end
    if nlim > 3
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,4,1),'-+k')
    end
    if nlim > 4
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,5,1),'-+y')
    end
    if nlim > 5
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,6,1),'-om')
    end
    if nlim > 6
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,7,1),'-^c')
    end

    subplot(1,2,2),loglog([1.0e2,1.0e0],[1.0e0,1.0e-2],'k-')
    subplot(1,2,2),loglog([1.0e2,1.0e0],[1.0e0,1.0e-4],'k-')
    grid off;
    xlabel('dx');
    title(title_str);
    axis(paxis);
    %axis square;
    %axis equal;
    hold off;

    % This section is just for the legend
    subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,1,1),'-+r')
    hold on;
    if nlim > 1
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,2,1),'-og')
    end
    if nlim > 2
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,3,1),'-^b')
    end
    if nlim > 3
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,4,1),'-+k')
    end
    if nlim > 4
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,5,1),'-+y')
    end
    if nlim > 5
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,6,1),'-om')
    end
    if nlim > 6
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,7,1),'-^c')
    end
    subplot(3,4,5),loglog([1.0e-1,1.0e-2],[1.0e-2,1.0e-3],'k--')
    subplot(3,4,5),loglog([1.0e-1,1.0e-2],[1.0e-2,1.0e-4],'k:')
    legend('None','Superbee','Lax-Wen','BeamWarm','Fromm','MinMod','MC','O(1)','O(2)','Location','northwest');
    grid on;
    axis([10 11 10 11]);
    axis off;
    hold off;
    print(out_file,'-color','-dpng');
  end
end

