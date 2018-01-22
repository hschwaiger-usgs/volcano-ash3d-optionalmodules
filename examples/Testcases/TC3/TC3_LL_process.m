clear all;

IsOct    = 0;
TestCase = 3;

% Set this to 1 if you want plots of the solution and errors.
% Set to 0 if you only want convergence plots and time plots.
cplotsol  = 0;
cploterr  = 1;
cploterr2 = 1;
cplotmass = 1;

nlim = 7;
nlim_label = char('LIM_NO','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_SB','LIM_MC');
nres  = char('1000','0500','0250');

resdx = [1.0000 0.500 0.250 0.125];
resnx = [72 144 288 566];
resny = [72 144 288 566];

xmin = 324.0;
xmax = 396.0;
ymin = -36.0;
ymax =  36.0;

nresmax  = 4;

GlobL1Errors  = zeros(nresmax,nlim);
GlobExecTimes = zeros(nresmax,nlim);
GlobalConserv = zeros(nresmax,nlim);
cinter = linspace(0.05,0.95,10);

figure(1);

  for il = 1:nlim
    for inres = 1:nresmax

      ifile = sprintf('DATA_LL/TC%i_LL__%s_%s_sol.dat',  ...
          TestCase,nlim_label(il,:),nres(inres,:));
      ierfile = sprintf('DATA_LL/TC%i_LL__%s_%s_err.dat',  ...
          TestCase,nlim_label(il,:),nres(inres,:));
      itimefile = sprintf('DATALL/TC%i_LL_%s_%s_etime.dat',  ...
          TestCase,nlim_label(il,:),nres(inres,:));

      dx = resdx(inres);
      nx = resnx(inres);
      ny = resny(inres);
      x_res = linspace(xmin,xmax,nx);
      y_res = linspace(ymin,ymax,ny);
      [x y] = meshgrid(x_res,y_res);

      % VARIABLES = "true","calc","err"
%      data = load(ifile);
%      truesol = reshape(data(:,1),nx,ny);
%      c       = reshape(data(:,2),nx,ny);
%      err     = reshape(data(:,3),nx,ny);
%      t_tot   = sum(sum(truesol));
%      c_tot   = sum(sum(c));
%      conserv = c_tot/t_tot;
      L1L2errors = load(ierfile);
%      truepeak = max(max(truesol));

      GlobL1Errors(inres,il) = L1L2errors(1);
      GlobL2Errors(inres,il) = L1L2errors(2);
      %GlobalConserv(inres,il)= conserv;
      GlobalConserv(inres,il)= L1L2errors(3);

      if (cplotsol == 1)
        if (inres == 1)
          figure(2,'Position',[10 10 1010 810])
        elseif (inres == nresmax)
          figure(1);
          subplot(3,4,il),contour(x,y,c,cinter);axis xy;axis([xmin xmax ymin ymax]);axis nolabel;
          figure(2);
        end
        limlab = regexprep(strtrim(nlim(il,:)), '_', '=');
        limlab = regexprep(strtrim(nlim(il,:)), '_', '=');
        %label = sprintf('%s\n%s\ndx = %f\nTotal L_1 Err = %f\nTotal L_2 Err = %f\nMass Cons Fac = %f', ...
        %   strtrim(nmeth(inmeth,:)),limlab,dx,L1L2errors(1),L1L2errors(2),conserv);
        end
        xlab = -0.5*(xmax-xmin)+xmin;
        ylab =  0.5*(ymax-ymin)+ymin;
        if (inres<4)
          subplot(4,4,(inres-1)*4+1),text(xlab,ylab,label);axis([xmin xmax ymin ymax]);axis off;
          subplot(4,4,(inres-1)*4+2),imagesc(x,y,truesol);axis xy;axis([xmin xmax ymin ymax]);title('True');
          subplot(4,4,(inres-1)*4+3),imagesc(x,y,c);axis xy;axis([xmin xmax ymin ymax]);title('Calc');
          subplot(4,4,(inres-1)*4+4),imagesc(x,y,abs(err));axis xy;axis([xmin xmax ymin ymax]);title('Abs(Error)');
        end
        if (inres == nresmax)
          ofile = sprintf('Solution_LL_%s.png',strtrim(nlim(il,:)));
          ofile2 = sprintf('Solution_LL_%s.eps',strtrim(nlim(il,:)));
          print(ofile,'-color','-dpng');
          print(ofile2,'-color','-depsc');
          close(2);
          if (il==nlim)
            figure(1);
            print('LL_compare.eps','-color','-depsc');
          end
        end
      end

    end

    if cplottime==1
      exectimes_raw = load(itimefile);
      GlobExecTimes(:,il) = exectimes_raw(:);
    end


close(1);

if cploterr==1
  figure;
  subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,1,1),'r-@')
  hold on;
  subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,2,1),'g-@')
  subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,3,1),'b-@')
  subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,4,1),'k-@')
  if nmethmax>1
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,1,2),'r-@o')
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,2,2),'g-@o')
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,3,2),'b-@o')
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,4,2),'k-@o')
  end
  if nmethmax>2
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,1,3),'r-@*')
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,2,3),'g-@*')
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,3,3),'b-@*')
    subplot(3,4,5),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,4,3),'k-@*')
  end
  grid off;
  axis([10 11 10 11]);
  axis off;
  legend('DCU-No','DCU-MM','DCU-SB','DCU-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL-1','SL-5','SL-10','SL-20',2)


  subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,1,1),'r-@')
  hold on;
  subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,2,1),'g-@')
  subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,3,1),'b-@')
  subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,4,1),'k-@')
  if nmethmax>1
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,1,2),'r-@o')
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,2,2),'g-@o')
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,3,2),'b-@o')
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,4,2),'k-@o')
  end
  if nmethmax>2
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,1,3),'r-@*')
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,2,3),'g-@*')
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,3,3),'b-@*')
    subplot(1,4,2),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,4,3),'k-@*')
  end
  subplot(1,4,2),loglog([1.0e1,1.0e-1],[1.0e-1,1.0e-3],'k-')
  subplot(1,4,2),loglog([1.0e1,1.0e-1],[1.0e-1,1.0e-5],'k-')

  grid on;
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC',1)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL2','SL4','SL8',4)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL5','SL10','SL20',4)
  xlabel('dx');
  title('L_1 Error');
  axis([0.1 10 1.0e-6 1.0e-1]);
  axis equal;
  hold off;

  %subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,1,1),'r-@')
  %hold on;
  %subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,2,1),'g-@')
  %subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,3,1),'b-@')
  %subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,4,1),'k-@')
  %if nmethmax>1
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,1,2),'r-@o')
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,2,2),'g-@o')
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,3,2),'b-@o')
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,4,2),'k-@o')
  %end
  %if nmethmax>2
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,1,3),'r-@*')
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,2,3),'g-@*')
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,3,3),'b-@*')
  %  subplot(1,2,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,4,3),'k-@*')
  %end
  %grid on;
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC',1)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL2','SL4','SL8',1)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL5','SL10','SL20',4)
  %xlabel('dx');
  %ylabel('L_2 Error');
  %axis equal;
  %hold off;

  subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,1,1),'r-@')
  hold on;
  subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,2,1),'g-@')
  subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,3,1),'b-@')
  subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,4,1),'k-@')
  if nmethmax>1
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,1,2),'r-@o')
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,2,2),'g-@o')
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,3,2),'b-@o')
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,4,2),'k-@o')
  end
  if nmethmax>2
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,1,3),'r-@*')
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,2,3),'g-@*')
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,3,3),'b-@*')
    subplot(1,4,3),loglog(resdx(1:nresmax),GlobalConserv(1:nresmax,4,3),'k-@*')
  end
  grid on;
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC',1)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL2','SL4','SL8',1)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL5','SL10','SL20',4)
  xlabel('dx');
  title('Mass Cons. Ratio');
  axis equal;
  %axis([0.1 10 1.0e-6 1.0e-1]);
  hold off;

  %print('Solution_Error_L1.png','-color','-dpng')
  %print('Solution_Error_L1.eps','-color','-depsc')
  %close();
end


if (cplottime == 1)
  %figure;
  %hold on;
  subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,1,1),'r-@')
  hold on;
  subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,2,1),'g-@')
  subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,3,1),'b-@')
  subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,4,1),'k-@')
  if nmethmax>1
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,1,2),'r-@o')
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,2,2),'g-@o')
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,3,2),'b-@o')
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,4,2),'k-@o')
  end
  if nmethmax>2
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,1,3),'r-@*')
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,2,3),'g-@*')
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,3,3),'b-@*')
   subplot(1,4,4),loglog(resdx(1:nresmax),GlobExecTimes(1:nresmax,4,3),'k-@*')
  end
  grid on;
  axis([0.1 10 1.0e-2 1.0e3]);
  axis equal;

  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC',1)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL2','SL4','SL8',1)
  %legend('DS-No','DS-MM','DS-SB','DS-MC','CTU-No','CTU-MM','CTU-SB','CTU-MC','SL1','SL5','SL10','SL20',1)
  xlabel('dx');
  title('Execution time (s)');
  hold off;
  print('Solution_ExecTime.png','-color','-dpng')
  print('Solution_ExecTime.eps','-color','-depsc')
  close()
end
