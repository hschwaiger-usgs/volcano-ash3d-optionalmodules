clear all;

IsOct    = 0;
TestCase = 3;

% Set this to 1 if you want plots of the solution and errors.
% Set to 0 if you only want convergence plots and time plots.
cplotsol  = 1;
cploterr1 = 1;
cploterr2 = 1;
cplotmass = 1;

nlim = 1;
nlim_label = char('LIM_NO','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_SB','LIM_MC');
nres  = char('50000','25000','12500','06250','03125');

resdx = [0.05000 0.025000 0.012500 0.006250 0.003125];
resnx = [40 80 160 320 640];
resny = [40 80 160 320 640];

xmin = -1.0;
xmax =  1.0;
ymin = -1.0;
ymax =  1.0;

nresmax  = 3;

GlobL1Errors  = zeros(nresmax,nlim);
GlobExecTimes = zeros(nresmax,nlim);
GlobalConserv = zeros(nresmax,nlim);
cinter = linspace(0.05,0.95,10);

figure(1);

  for il = 1:nlim
    for inres = 1:nresmax

      ifile = sprintf('DATA_xy/TC%i_XY_%s_%s_sol.dat',  ...
          TestCase,nlim_label(il,:),nres(inres,:));
      ierfile = sprintf('DATA_xy/TC%i_XY_%s_%s_err.dat',  ...
          TestCase,nlim_label(il,:),nres(inres,:));
      itimefile = sprintf('DATA_xy/TC%i_XY_%s_etime.dat',  ...
          TestCase,nlim_label(il,:));

      dx = resdx(inres);
      nx = resnx(inres);
      ny = resny(inres);

      x_res = linspace(-1.0+0.5*dx,1.0-0.5*dx,nx);
      y_res = linspace(-1.0+0.5*dx,1.0-0.5*dx,ny);
      [xg yg] = meshgrid(x_res,y_res);
      x =  xg(1,:);
      y =  yg(:,1);

      % VARIABLES = "true","calc","err"
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
      GlobalConserv(inres,il)= L1L2errors(3);

      if (cplotsol == 1)
        if (inres == 1)
          %figure(2,'Position',[10 10 1010 810])
          figure(2);
        elseif (inres == nresmax)
          figure(1);
          %subplot(3,4,inlim),imagesc(x,y,c);axis xy;axis([xmin xmax ymin ymax]);axis off;
          subplot(3,4,il),contour(x,y,c,cinter);axis xy;axis([xmin xmax ymin ymax]);%axis nolabel;
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

        if (inres == nresmax)
          ofile = sprintf('PLOTS/Solution_%s.png',  ...
            nlim_label(il,:));
          ofile2 = sprintf('PLOTS/Solution_%s.eps',  ...
            nlim_label(il,:));
          %if IsOct == 1
          %  %print(ofile,'-color','-dpng');
            print(ofile2,'-color','-depsc');
          %else
          %  print(2,'-depsc2',ofile2)
          %end
          close(2);
          if (il==nlim)
            figure(1);
            %if IsOct == 1
<<<<<<< HEAD
              print('PLOTS/xyz_compare.eps','-color','-depsc');
=======
              print('xyz_compare.eps','-color','-depsc');
>>>>>>> f628dc1caa97b0b7f2f4f7e9241a5b90f92a6629
            %else
            %  print(1,'-depsc2','xyz_compare.eps')
            %end
          end
        end
      end

    end

  end

close(1);

plotflags = [cploterr1, cploterr2, cplotmass];
for iplot = 1:3
  if plotflags(iplot)==1
    if iplot==1
      plotvar = GlobL1Errors;
      title_str = 'L_1 Error';
      out_file = 'Solution_Error_L1.eps';
      paxis = [1.0e-3 1.0e-1 1.0e-3 1.0e-1];
    elseif iplot==2
      plotvar = GlobL2Errors;
      title_str = 'L_2 Error';
      out_file = 'Solution_Error_L2.eps';
      paxis = [1.0e-3 1.0e-1 1.0e-2 1.0e0];
    elseif iplot==3
      plotvar = GlobalConserv;
      title_str = 'Mass Cons. Ratio';
      out_file = 'Solution_Error_MassCons.eps';
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

    subplot(1,2,2),loglog([1.0e-1,1.0e-3],[1.0e-1,1.0e-3],'k-')
    subplot(1,2,2),loglog([1.0e-1,1.0e-3],[1.0e-1,1.0e-5],'k-')
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
    subplot(3,4,5),loglog([1.0e-1,1.0e-2],[1.0e-2,1.0e-3],'k--')
    subplot(3,4,5),loglog([1.0e-1,1.0e-2],[1.0e-2,1.0e-4],'k:')
    legend('None','Lax-Wen','BeamWarm','Fromm','MinMod','Superbee','MC','O(1)','O(2)','Location','northwest');
    grid on;
    axis([10 11 10 11]);
    axis off;
    hold off;



    %if IsOct == 1
      print(out_file,'-color','-depsc');
      %close();
    %else
    %  print(1,'-depsc2',out_file);
    %  close();
    %end
  end
end

