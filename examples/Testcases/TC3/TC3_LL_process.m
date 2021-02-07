clear all;

TestCase = 3;

% Set this to 1 if you want plots of the solution and errors.
% Set to 0 if you only want convergence plots and time plots.
cplotsol  = 1;
cploterr  = 0;

nlim = 7;
nlim_label = char('LIM_NO','LIM_LW','LIM_BW','LIM_FM','LIM_MM','LIM_SB','LIM_MC');
nres  = char('2000','1000','0500','0250','0125');

resdx = [2.0000 1.0000 0.500 0.250 0.125];
resnx = [36 72 144 288 576];
resny = [36 72 144 288 576];

xmin = 324.0;
xmax = 396.0;
ymin = -36.0;
ymax =  36.0;

nresmax  = 5;

GlobL1Errors  = zeros(nresmax,nlim);
GlobExecTimes = zeros(nresmax,nlim);
GlobalConserv = zeros(nresmax,nlim);
cinter = linspace(0.05,0.95,10);

if cplotsol==1
  il = 7;
  figure(1,'Position',[10 10 1010 210])
  for inres = 1:nresmax

    ifile = sprintf('DATA_LL/TC%i_LL__%s_%s_sol.dat',  ...
        TestCase,nlim_label(il,:),nres(inres,:));

    dx = resdx(inres);
    nx = resnx(inres);
    ny = resny(inres);
    x_res = linspace(xmin,xmax,nx);
    y_res = linspace(ymin,ymax,ny);
    [x y] = meshgrid(x_res,y_res);

    % VARIABLES = "true","calc","err"
    data = load(ifile);
    truesol = reshape(data(:,1),nx,ny);
    c       = reshape(data(:,2),nx,ny);
    err     = reshape(data(:,3),nx,ny);
    t_tot   = sum(sum(truesol));
    c_tot   = sum(sum(c));
    conserv = c_tot/t_tot;

    subplot(2,nresmax+1,inres),contour(x,y,c,cinter);axis xy;
          axis([xmin xmax ymin ymax]);axis square ;axis nolabel;
    if (inres == nresmax)
      subplot(2,nresmax+1,nresmax+1),imagesc(x,y,truesol);axis xy;
          axis([xmin xmax ymin ymax]);axis square ;axis nolabel;title('True');
    end
    subplot(2,nresmax+1,inres+nresmax+1),contour(x,y,abs(err));axis xy;
          axis([xmin xmax ymin ymax]);axis square ;axis nolabel;
  end
  print('TC3_LL_sol_LIM_.eps','-color','-depsc');
end

if cploterr==1
  for il = 1:nlim
    for inres = 1:nresmax
      ierfile = sprintf('DATA_LL/TC%i_LL__%s_%s_err.dat',  ...
        TestCase,nlim_label(il,:),nres(inres,:));
      L1L2errors = load(ierfile);
      GlobL1Errors(inres,il) = L1L2errors(1);
      GlobL2Errors(inres,il) = L1L2errors(2);
      GlobMCErrors(inres,il)= L1L2errors(3);
    end
  end
  figure;
  %Specify point of O(1) and O(2) intersection
  Ox=1.0e0;
  Oy=1.0e-1;
  dxmin=1.0e-1;
  dxmax=1.0e0;
  L1min=1.0e-3;
  L1max=1.0e-1;
  L2min=1.0e-6;
  L2max=1.0e-4;
  Ermin=1.0e-16;
  Ermax=1.0e-13;
  Ermin=1.0e-16;
  Ermax=1.0e0;
  
  subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-2],'k--');hold on;
  subplot(2,3,1),loglog([Ox,Ox*1.0e-2],[Oy,Oy*1.0e-4],'k:')
  subplot(2,3,1),loglog(resdx(1:nresmax),GlobL1Errors(1:nresmax,:));axis([dxmin dxmax L1min L1max]);hold
  ylabel('L_1 Error')
  subplot(2,3,2),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,:));axis([dxmin dxmax L2min L2max]);
  ylabel('L_2 Error')
  xlabel('dx')
  subplot(2,3,3),loglog(resdx(1:nresmax),GlobMCErrors(1:nresmax,:));axis([dxmin dxmax Ermin Ermax]);
  ylabel('MC Error')
  xlabel('dx')

  subplot(2,3,5),loglog(resdx(1:nresmax),GlobL2Errors(1:nresmax,:));axis([dxmin dxmax 1.0e-8 1.0e-6]);
  legend('NO','LW','BW','FM','MM','SB','MC')

  print('TC3_LL_convergence.eps','-color','-depsc')
  %close();
end
