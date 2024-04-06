function Simulations_1DOG(in_K, in_model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  "Mechanical models of pattern and form in biological tissues:      %%%
%%%       the role of stress-strain constitutive equations"             %%%
%%%                                                                     %%%
%%%      C. Villa (*), M.A.J. Chaplain, A. Gerisch (**), T. Lorenzi     %%%
%%%                                                                     %%%
%%%            Bullettin of Mathematical Biology (2021)                 %%%
%%%                                                                     %%%
%%%                                                                     %%%
%%% (*) cv23[at]st-andrews.ac.uk                                        %%%
%%% (**) gerisch[at]mathematik.tu-darmstadt.de                          %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  1D Simulations (Kelvin-Voigt and Maxwell models)                   %%%
%%%  For details about the equations and schemes please see the         %%%
%%%  manuscript indicated above and the supplementary material          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function Simulations_1D(in_K, in_model)
% Runs the spatially one-dimensional simulation using viscoelastic 
% model in_model (either 'Maxwell' or 'Kelvin-Voigt') on a uniform grid
% with in_K grid cells. 
% The results are saved in file saved_y1D_[in_model]_[in_K].mat and 
% a video in saved_y1D_[in_model]_[in_K].avi .
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulations_1D.m: simulates a 1D mechanical model of pattern formation%
%%% Copyright (C) 2021 C. Villa and A. Gerisch                            %
%%%                                                                       %
%%% This program is free software: you can redistribute it and/or modify  %
%%% it under the terms of the GNU General Public License as published by  %
%%% the Free Software Foundation, either version 3 of the License, or     %
%%% (at your option) any later version.                                   %
%%%                                                                       %
%%% This program is distributed in the hope that it will be useful,       %
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%%% GNU General Public License for more details.                          %
%%%                                                                       %
%%% You should have received a copy of the GNU General Public License     %
%%% along with this program.  If not, see <https://www.gnu.org/licenses/>.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
close all

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')


%% Numerical set up
par.seed = 44;                 % Seed for random number generator
par.model = in_model;          % type of viscoelastic model
%par.model = 'Maxwell';
%par.model = 'Kelvin-Voigt'; 
par.K = in_K;                  % Number of spatial grid cells
par.L = 1;                     % Domain length
x = linspace(0,par.L,par.K+1); % Discretise spatial domain
par.dx = x(2)-x(1);            % Cell size
t0 = 0;                        % Initial time
tf = 1.5;                    % Final time
tspan = linspace(t0,tf,200);   % Time span

%% Initial conditions - eq.(28)
steadystate = [ones(2*par.K,1); zeros(par.K,1)];
f = 1+exp(-(x(1:end-1)-0.5*par.L).^2);        % adding a gaussian noise
n0 = (f.*ones(1,par.K)).';      % cells
p0 = (f.*ones(1,par.K)).';      % collagen
u0 = zeros(par.K,1);            % displacement
y0 = [n0;p0;u0];                % concactanate arrays

%% Solve with ODE15i
%%% Compute consistent yp0 
%res = @(y,yp)(norm(mechanochemical(y,yp,par),'inf'));
res = @(y,yp)(norm(mechanochemical(y,yp,par)));
disp(['residuum of steady state = ' ...
    num2str(res(steadystate,0*steadystate), '%15.10e')]);
[y0,yp0,resnorm] = decic(@(t,y,yp)(mechanochemical(y,yp,par)), t0, ...
    y0, zeros(3*par.K,1), zeros(3*par.K,1), zeros(3*par.K,1));
disp(['residuum (from decic) of IC = ' num2str(resnorm, '%15.10e')]);
disp(['residuum (from res()) of IC = ' num2str(res(y0,yp0), '%15.10e')]);
%%% Solve
tic
% options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',1e-6);
[t,y] = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
toc
%%% Save computed solution to file
filename = ['saved_y1D_' par.model '_' num2str(par.K)];
save(filename, 't', 'y', 'par', 'x');

%% Plot
video_on = true; % Record video: YES (true), NO (false)
video_filename = [filename '.avi'];
plot_solution(x,y,t,par,video_on,video_filename);

end

%% Main function implementing the model
function f = mechanochemical(y,yp,par)

    %% Parameter values
    eta = 1;      % viscosity
    E = 1;        % elasticity
    D = 0.01;     % diffusion
    alpha = 0;    % haptotaxis
    r = 0;        % proliferation
    an = 0.5;     % cell growth
    s = 0.1;      % substrate elasticity
    beta = 0;     % long range traction 
    lambda = 0;   % (cell traction) saturation coefficient
    tau = 0.2;    % cell traction 
    dn = 0.5;     % cell death rate
    m = 0.5;      % collagen production rate
    dp = 0.5;     % collagen decay rate
    k = 5;      % stress induced by cells hills function parameter

    %%% Choose constitutive model (see Table 1)
    switch par.model
        case 'Kelvin-Voigt' % Kelvin Voigt  - eq.(3)
            [a0,a1,b0,b1] = deal(1,0,E,eta);
        case 'Maxwell' % Maxwell - eq.(4)
            [a0,a1,b0,b1] = deal(1/eta,1/E,0,1);
        otherwise
            error('Unknown constitutive model')
    end
    
    %%% Reshape input vectors
    [n,p,u] = deal(y(1:par.K),y(par.K+1:2*par.K),y(2*par.K+1:3*par.K));
    [np,pp,up] = deal(yp(1:par.K),yp(par.K+1:2*par.K),...
        yp(2*par.K+1:3*par.K));
    
    %%% Equation for n
    % Advection velocity at grid cell interfaces - eq.(S.6)
    up_Avx1 = Avx1(up); 
    vx1 = alpha*Mx1_face(p,par) + up_Avx1; 
    % fn(n,n',p,u') = 0 - eq.(S.5)
    % fn = np - D*Mxx1(n,par) + MA1(n, vx1, par) - r*n.*(1-n);
    fn = np - D*Mxx1(n,par) + MA1(n, vx1, par) + an*ones(size(n)) - dn*n - r*n.*(1-n);

    %%% Equation for p
    % fp(p,p',u') = 0 - eq.(S.10)
    fp = pp + MA1(p, up_Avx1, par) + m*n./(1+n) - dp*p;
   
    %%% Equation for u 
    % Traction term - eq.(S.12)-(S.14)
    pexp = 2; 
    fn1 = n./(1+lambda*n.^pexp); % Lambda_1
    fn2 = ((1-(pexp-1)*lambda*n.^pexp)./((1+lambda*n.^pexp).^2)); %Lambda_2
    fn3 = (n.^k)./(1+n.^k); % hills function
    fp1 = p + beta*Mxx1(p, par); % M_T1 P
    fp2 = pp + beta*Mxx1(pp, par); % M_T1 P'
    Tr = tau*(a0*fn3.*fp1 + a1*(fn2.*np.*fp1+fn1.*fp2));
    % fu(n,n',p,p',u,u') = 0 - eq.(S.11)
    fu = b1*Mxx1(up, par) + b0*Mxx1(u, par) ...
        + Mx1_center(Tr, par) - a1*s*(p.*up + pp.*u) - a0*s*p.*u;
    
    %%% Full system - eq.(S.1)
    f = [fn; fp; fu];
end
    

%% Annexed functions

%%% Compute variable at grid cell interfaces - def.(S.3)
function avx1 = Avx1(y)
    avx1 = 0.5*(y +y([2:end,1]));
end

%%% Compute first order derivative on grid cell interfaces 
%%% (second order approximation, central) - def.(S.4)
function dx1 = Mx1_face(y, par)
    persistent Mdx;
    if size(Mdx,2) ~= size(y,1)
      disp('Mx1_face  : resetting persistent matrix.');  
      c = (1/par.dx)*sparse([-1; zeros(par.K-2,1); 1]); % First column
      Mdx = toeplitz(c,[c(1), c(end:-1:2)']); 
    end
    dx1 = Mdx*y;
end

%%% Compute first order derivative on grid cell centers
%%% (second order approximation, central) - def.(S.2)
function dx1 = Mx1_center(y, par)
    persistent Mdx;
    if size(Mdx,2) ~= size(y,1)
      disp('Mx1_center: resetting persistent matrix.');  
      c = (1/(2*par.dx))*sparse([0;-1; zeros(par.K-3,1); 1]); 
      Mdx = toeplitz(c,[c(1), c(end:-1:2)']); 
    end
    dx1 = Mdx*y;
end

%%% Compute second order derivative on grid cell centers
%%% (second order approximation, central) - def.(S.2)
function dxx1 = Mxx1(y, par)
    persistent Mdxx;
    if size(Mdxx,2) ~= size(y,1)
      disp('Mxx1      : resetting persistent matrix.');  
      c = (1/(par.dx^2))*sparse([-2; 1; zeros(par.K-3,1); 1]);
      Mdxx = toeplitz(c,[c(1), c(end:-1:2)']);
    end
    dxx1 = Mdxx*y;    
end

%%% Compute advection at grid cell centers using first order upwinding
%%% with advective velocity given at grid cell interfaces - def.(S.7)
function fluxdiffx1 = MA1(y, vel, par)
    % find pos/neg velocity entries
    velPos_ind = (vel>0);
    velNeg_ind = (~velPos_ind);
    % compute flux accross cell interfaces using first order upwinding 
    % def.(S.8)-(S.9)
    flux = NaN(size(y));
    flux(velPos_ind) = y(velPos_ind).*vel(velPos_ind);
    yshift = y([2:end 1]);
    flux(velNeg_ind) = yshift(velNeg_ind).*vel(velNeg_ind);
    % add entries for flux along x1=0 by periodicity
    flux = flux([end,1:end]);
    % compute flux difference per grid cell - def.(S.7) and (S.4)
    fluxdiffx1 = (1/par.dx)*(flux(2:end)-flux(1:(end-1)));
end

%%% Plot solution 
function plot_solution(x,y,t,par,video_on,video_filename)
    if video_on % Initialise video
        vid = VideoWriter(video_filename);
        open(vid);
        figure('Units','normalized','Position',[0 0 0.5 0.45])
    end
    maxu = 10^(-6);
    for i=1:length(t)
        clf
        n = [y(i,1*par.K),y(i,1:par.K)];
        p = [y(i,2*par.K),y(i,par.K+1:2*par.K)];
        u = [y(i,3*par.K),y(i,2*par.K+1:3*par.K)];
        if max(abs(u))>maxu
          maxu = max(abs(u))+0.1*max(abs(u));
        end
        subplot(1,3,1)
        plot(x,n)
        title('$n(t,x)$')
        ylim([0,max(2,max(n))])
        axis square
        subplot(1,3,2)
        plot(x,p)
        title('$\rho(t,x)$')
        ylim([0,max(2,max(p))])
        axis square
        subplot(1,3,3)
        plot(x,u)
        title('$u(t,x)$')
        axis square
        ylim([-maxu,maxu])
        a = axes;
        t1 = title([par.model ' (t=',num2str(t(i)),')']);
        a.Visible = 'off'; 
        t1.Visible = 'on'; 
        drawnow
        if video_on % Record video
            frame = getframe(gcf);
            size(frame.cdata);
            writeVideo(vid,frame);
            pause(0.1)
        end
    end
    if video_on % Close video
        close(vid)
    end
end
