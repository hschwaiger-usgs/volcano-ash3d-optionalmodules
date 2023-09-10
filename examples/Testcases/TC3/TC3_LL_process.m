#!/usr/bin/octave -qf
clear all;

IsOct    = 0;
TestCase = 3;

% Load parameters from run script
nresmax = load('TC3_LL_idx.dat');
datlim  = load('TC3_LL_lim.dat');
limmin=datlim(1)+1;
limmax=datlim(2)+1;

% Set this to 1 if you want plots of the solution and errors.
% Set to 0 if you only want convergence plots and time plots.
cplotsol  = 1; % Plot solution (row=true, calc, err, col=dx)
cploterr1 = 1; % Plot L1 error
cploterr2 = 0; % Plot L2 error
cplotmass = 1; % Plot Mass Consv error

resdx = [2.0000 1.0000 0.500 0.250 0.125];
resnx = [36 72 144 288 576];
resny = [36 72 144 288 576];

nlim = limmax-limmin+1;
nlim_label = char('LIM_NO','LIM_SB','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_MC');
nlimLegLabel=char('NO','SB','LW','BW','FM','MM','MC');
nres  = char('2000','1000','0500','0250','0125');

xmin = 324.0;
xmax = 396.0;
ymin = -36.0;
ymax =  36.0;

GlobL1Errors  = zeros(nresmax,nlim);
GlobL2Errors  = zeros(nresmax,nlim);
GlobMCErrors  = zeros(nresmax,nlim);
GlobExecTimes = zeros(nresmax,nlim);

cinter = linspace(0.05,0.95,10);

figure(1);  % Figure 1 will be for the final contour plot of the individual runs
for il = 1:nlim
  for inres = 1:nresmax

    ifile = sprintf('DATA/TC%i_LL_%s_%s_sol.dat',  ...
        TestCase,nlim_label(il,:),nres(inres,:));
    ierfile = sprintf('DATA/TC%i_LL_%s_%s_err.dat',  ...
        TestCase,nlim_label(il,:),nres(inres,:));
    itimefile = sprintf('DATA/TC%i_LL_%s_etime.dat',  ...
        TestCase,nlim_label(il,:));

    dx = resdx(inres);
    nx = resnx(inres);
    ny = resny(inres);

    x_res = linspace(xmin,xmax,nx);
    y_res = linspace(ymin,ymax,ny);
    [x y] = meshgrid(x_res,y_res);

    % Reading in solution file with columns of "true","calc","err"
    data = load(ifile);
    truesol = reshape(data(:,1),nx,ny);
    c       = reshape(data(:,2),nx,ny);
    err     = reshape(data(:,3),nx,ny);
    t_tot   = sum(sum(truesol));
    c_tot   = sum(sum(c));
    conserv = c_tot/t_tot;
    L1L2errors = load(ierfile);
    truepeak = max(max(truesol));

    GlobL1Errors(inres,il) = L1L2errors(1);
    GlobL2Errors(inres,il) = L1L2errors(2);
    GlobMCErrors(inres,il) = L1L2errors(3);

    if (cplotsol == 1)
      if (inres == 1)
        figure(2);  % Figure 2 will be for the solution plots, if requested
      elseif (inres == nresmax)
        figure(1);
        subplot(3,3,il),contour(x,y,c,cinter);axis xy;axis([xmin xmax ymin ymax]);axis nolabel;
        title(nlimLegLabel(il,1:2));
        figure(2);
      end
      limlab = regexprep(nlim_label(il,:), '_', '=');
      limlab = regexprep(nlim_label(il,:), '_', '=');
      label = sprintf('%s\n%s\ndx = %f\nTotal L_1 Err = %f\nTotal L_2 Err = %f\nMass Cons Fac = %f', ...
         limlab,dx,L1L2errors(1),L1L2errors(2),conserv);

      subplot(nresmax,4,(inres-1)*4+1),text(-0.5,0.5,label);axis([xmin xmax ymin ymax]);axis off;
      subplot(nresmax,4,(inres-1)*4+2),imagesc(x,y,truesol);axis xy;axis([xmin xmax ymin ymax]);title('True');
      subplot(nresmax,4,(inres-1)*4+3),imagesc(x,y,c);axis xy;axis([xmin xmax ymin ymax]);title('Calc');
      subplot(nresmax,4,(inres-1)*4+4),imagesc(x,y,abs(err));axis xy;axis([xmin xmax ymin ymax]);title('Abs(Error)');
       % If this is the last dx, plot the solution
      if (inres == nresmax)
        ofile = sprintf('PLOTS/TC3_LL_Solution_LL_%s.png',  ...
          nlim_label(il,:));
        print(ofile,'-color','-dpng');
        if (il==nlim)
          figure(1);
          print('PLOTS/TC3_LL_compare.png','-color','-dpng');
        end
        close(2);
      end
    end

  end
end
close(1);
% Write out the L1 errors for each limiter
ConvRate=zeros(nlim,1);
ConvRate(1:nlim)=log(GlobL1Errors(1,1:nlim)./GlobL1Errors(nresmax,1:nlim))/log(resdx(1)/resdx(nresmax));

save ("-ascii","DATA/TC3_ConvRate_LL.dat","ConvRate")

plotflags = [cploterr1, cploterr2, cplotmass];
for iplot = 1:3
  if plotflags(iplot)==1
    if iplot==1
      plotvar = GlobL1Errors;
      title_str = 'L_1 Error';
      out_file = 'PLOTS/TC3_LL_Solution_Error_L1.png';
      paxis = [1.0e-3 1.0e-1 1.0e-3 1.0e-1];
    elseif iplot==2
      plotvar = GlobL2Errors;
      title_str = 'L_2 Error';
      out_file = 'PLOTS/TC3_LL_Solution_Error_L2.png';
      paxis = [1.0e-3 1.0e-1 1.0e-2 1.0e0];
    elseif iplot==3
      plotvar = GlobMCErrors;
      title_str = 'Mass Cons. Ratio';
      out_file = 'PLOTS/TC3_LL_Solution_Error_MassCons.png';
      paxis = [1.0e-3 1.0e-1 1.0e-7 1.0e-1];
    end
    figure;
    subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,1,1),'-+r')
    hold on;
    if nlim > 1
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,2,1),'-+g')
    end
    if nlim > 2
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,3,1),'-+b')
    end
    if nlim > 3
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,4,1),'-+k')
    end
    if nlim > 4
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,5,1),'-+y')
    end
    if nlim > 5
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,6,1),'-+m')
    end
    if nlim > 6
      subplot(1,2,2),loglog(resdx(1:nresmax),plotvar(1:nresmax,7,1),'-+c')
    end

    subplot(1,2,2),loglog([1.0e1,1.0e0],[1.0e-1,1.0e-2],'k-')
    subplot(1,2,2),loglog([1.0e1,1.0e0],[1.0e-1,1.0e-3],'k-')

    grid off;
    xlabel('dx');
    title(title_str);
    %axis(paxis);
    %axis square;
    %axis equal;
    hold off;

    % This section is just for the legend
    subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,1,1),'-+r')
    hold on;
    if nlim > 1
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,2,1),'-+g')
    end
    if nlim > 2
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,3,1),'-+b')
    end
    if nlim > 3
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,4,1),'-+k')
    end
    if nlim > 4
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,5,1),'-+y')
    end
    if nlim > 5
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,6,1),'-+m')
    end
    if nlim > 6
      subplot(3,4,5),loglog(resdx(1:nresmax),plotvar(1:nresmax,7,1),'-+c')
    end

    subplot(3,4,5),loglog([1.0e1,1.0e-1],[1.0e0,1.0e-2],'k--')
    subplot(3,4,5),loglog([1.0e1,1.0e-1],[1.0e0,1.0e-3],'k:')

    legend('None','Superbee','Lax-Wen','BeamWarm','Fromm','MinMod','MC','O(1)','O(2)','Location','northwest');
    grid on;
    axis([10 11 10 11]);
    axis off;
    hold off;
    print(out_file,'-color','-dpng');
  end
end

