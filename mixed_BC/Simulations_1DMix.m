function Simulations_1DMix(in_K)
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
par.K = in_K;                    % Number of spatial grid cells
par.L = 1;                       % Domain length
par.Nright = 1;
par.Pright = 1;
par.Uleft = 0;
par.Uright = 0;
x = linspace(0,par.L,par.K);
par.x = linspace(0,par.L,par.K); % Discretise spatial domain
par.dx = par.x(2)-par.x(1);      % Cell size
t0 = 0;                          % Initial time
tf = 200;                         % Final time
tspan = linspace(t0,tf,100);     % Time span
%% Initial conditions - eq.(28)
steadystate = [ones(2*par.K,1); zeros(par.K,1)];
f = 1+0.5*exp(-((par.x)./0.2).^2);
% f1 = 2 + 10*sin(0.1*pi.*par.x);
% f2 = 1 + 0.5*par.x;
% f3 = 2 + 10*cos(0.05*pi.*par.x);
n0 = (f.*ones(1,par.K)).';      % cells
% n0 = ones(par.K,1);
% p0 = (f.*ones(1,par.K)).';      % collagen
p0 = ones(par.K,1);
u0 = zeros(par.K,1);            % displacement
% u0 = sin(pi.*x+pi/2);
% u0 = (f2.*ones(1,par.K)).';
y0 = [n0;p0;u0];

%% Solve with ODE15i
%%% Compute consistent yp0 
%res = @(y,yp)(norm(mechanochemical(y,yp,par),'inf'));
res = @(y,yp)(norm(mechanochemical(y,yp,par)));
disp(['residuum of steady state = ' ...
    num2str(res(steadystate,0*steadystate), '%15.10e')]);
opt = odeset('RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));
nfixed = zeros(par.K,1);
pfixed = zeros(par.K,1);
ufixed = zeros(par.K,1);
[y0,yp0,resnorm] = decic(@(t,y,yp)(mechanochemical(y,yp,par)), t0, ...
    y0, [nfixed; pfixed; ufixed], ...
    [zeros(3*par.K,1)], zeros(3*par.K,1), opt);
disp(['residuum (from decic) of IC = ' num2str(resnorm, '%15.10e')]);
disp(['residuum (from res()) of IC = ' num2str(res(y0,yp0), '%15.10e')]);
%% Solve
tic
[t,y] = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
toc
%%% Save computed solution to file
filename = ['saved_y1D_' num2str(par.K)];
save(filename, 't', 'y', 'par', 'x');

%% Plot
video_on = true; % Record video: YES (true), NO (false)
video_filename = [filename '.avi'];
plot_solution(x,y,t,par,video_on,video_filename);
end

%% Main function implementing the model
function f = mechanochemical(y,yp,par)
    %%% Parameter values
    eta = 1;      % viscosity
    E = 1;        % elasticity
    D = 0.01;     % diffusion
    r = 0;        % proliferation
    s = 1;        % substrate elasticity
    tau = 0.5;   % cell traction 
    an = 1;       % cell recruitment rate
    dn = 1;       % cell decay rate
    m = 1;      % collagen production rate
    dp = 1;     % collagen decay rate
    a1 = 0;     % stress related recruitment rate
    
    %%% Reshape input vectors
    [n,p,u] = deal(y(1:par.K),y(par.K+1:2*par.K),y(2*par.K+1:3*par.K));
    ntilde = [n(2); n; par.Nright]; 
    ptilde = [p(2); p; par.Pright]; 
    utilde = [par.Uleft; u; par.Uright];
    [np,pp,up] = deal(yp(1:par.K),yp(par.K+1:2*par.K),...
         yp(2*par.K+1:3*par.K));
    uptilde = [0; up; 0];
    
    %%% Equation for n
    % Advection velocity at grid cell interfaces - eq.(S.6)
    Tr = tau*p.*n;
    sig = eta*Mx(uptilde, par) + E*Mx(utilde,par) + Tr;
    % fn(n,n',p,u') = 0 - eq.(S.5)
    fn = np - D*Mxx(ntilde, par) + MA1(ntilde, up, par) + a1*sig - an*ones(size(n)) + dn*n - r*n.*(1-n);

    %%% Equation for p
    % fp(p,p',u') = 0 - eq.(S.10)
    fp = pp + MA1(ptilde, up, par) - m*n + dp*p;

    %%% Equation for u 
    % Traction term - eq.(S.12)-(S.14)
    n0 = 1.5;
    k1 = 10;
    h1 = (n.^k1)./(n0.^k1 + n.^k1);
    Tr = tau*p.*n;
    Trtilde = [Tr(2); Tr; tau*ptilde(end)*ntilde(end)];
    % Trtilde = [tau*p(2)*n(2); Tr; tau*ptilde(end)*ntilde(end)];
    % Trtilde = [Tr(2); tau*conv(p,n, "same"); tau*ptilde(end)*ntilde(end)];
    % I don't know why but when I have this scheme for traction force the
    % displacment BCs no longer gets distorted. 
    % Trtilde = [Tr(2); Tr; Tr(end-1)];
    % fu(n,n',p,p',u,u') = 0 - eq.(S.11)
    fu = eta*Mxx(uptilde, par) + E*Mxx(utilde,par) + Trx(Trtilde, par, tau) - s*p.*u;
    % fu = up - 0.1*sin(pi.*par.x+pi/2).';

    %%% Full system - eq.(S.1)
    f = [fn; fp; fu];
end

%% Annexed functions

%%% Compute first order derivative on grid cell interfaces 
%%% second order approximation, central
function dx1 = Mx(y,par)
    persistent Mdx;     % Mdx is a K x K+2 matrix
    c1 = [-1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdx = (0.5/par.dx)*[c1, diag(ones(par.K-1,1),1)-diag(ones(par.K-1,1),-1) ,cn];
    dx1 = Mdx*y;
end

%%% Compute second order derivative on grid cell interfaces 
%%% second order approximation, central
function dx2 = Mxx(y,par)
    persistent Mdxx;     % Mdx is a K x K+2 matrix
    c1 = [1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdxx = (1/par.dx/par.dx)*[c1, -2*eye(par.K) + diag(ones(par.K-1,1),1) + diag(ones(par.K-1,1),-1) ,cn];
    dx2 = Mdxx*y;
end

function trgradient = Trx(Trtilde, par, tau)
    trgradient = Mx(Trtilde, par);
    trgradient(1) = 0;
    trgradient(end) = (0.5/par.dx)*(tau-Trtilde(end-2));
end

%%% Compute advection at grid cell interfaces using first order upwinding
%%% with advective velocity given at grid cell interfaces - def.(S.7)
function fluxdiffx1 = MA1(y, vel, par)
    % compute flux accross cell interfaces using first order upwinding
    rightBC = y(end);
    y = y(2:end-1);
    % flux at K grid interfaces
    flux = NaN(size(vel));

    % take the cell center values for material, using the K by K-1 matrix
    % get K-1 grid cell centers
    c1 = [1; zeros(par.K-2,1)];
    yavg = 0.5*[c1, eye(par.K-1)+diag(ones(par.K-2,1),-1)]*y;
    for i = 2:size(flux)-1
        if (vel(i)>0)
            flux(i) = vel(i)*yavg(i-1);
        else 
            flux(i) = vel(i)*yavg(i);
        end
    end

    % for flux at the left Neumann boundary
    if vel(1)<0
    flux(1) = vel(1)*yavg(1);
    else
    flux(1) = vel(1)*y(1);
    end
    
    % for flux at the right boundary
    if vel(end)>0
        flux(end) = vel(end)*yavg(end);
    else
        flux(end) = vel(end)*rightBC;
    end

    % compute flux difference per grid cell - def.(S.7) and (S.4)
    fluxdiffx1 = Mx([flux(2); flux; flux(end-1)],par);
    fluxdiffx1(1) = (flux(2)-flux(1))/par.dx;
    fluxdiffx1(end) = (flux(end)-flux(end-1))/par.dx;
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
        n = [y(i,1:par.K)];
        p = [y(i,par.K+1:2*par.K)];
        u = [y(i,2*par.K+1:3*par.K)];
        if max(abs(u))>maxu
          maxu = max(abs(u))+0.1*max(abs(u));
        end
        subplot(1,3,1)
        plot(x,n)
        title('$n(t,x)$')
        % ylim([0,1.5])
        % ylim([0,max(2,max(n))])
        % ylim([0,12])
        axis square
        subplot(1,3,2)
        plot(x,p)
        title('$\rho(t,x)$')
        % ylim([0,1.5])
        % ylim([0,max(2,max(p))])
        axis square
        subplot(1,3,3)
        plot(x,u)
        title('$u(t,x)$')
        axis square
        % ylim([-maxu,maxu])
        % subplot(1,4,4)
        % plot(x,Sig)
        % title('$Sig(t,x)$')
        % axis square
        % ylim([min(Sig), max(Sig)])
        a = axes;
        t1 = title([' (t=',num2str(t(i)),')']);
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

%%% Plot transport rho, u, v
function plot_transport(x,y,t,par,video_on,video_filename)
    if video_on % Initialise video
        vid = VideoWriter(video_filename);
        open(vid);
        figure('Units','normalized','Position',[0 0 0.5 0.45])
    end
    for i=1:length(t)
        clf
        p = [y(i,par.K+1:2*par.K)];
        u = [y(i,2*par.K+1:3*par.K)];
        v = 0.1*sin(pi.*par.x+pi/2);
        subplot(1,3,1)
        plot(x,p)
        title('$p(t,x)$')
        axis square
        subplot(1,3,2)
        plot(x,u)
        title('$u(t,x)$')
        axis square
        subplot(1,3,3)
        plot(x,v)
        title('$v(x)$')
        axis square
        a = axes;
        t1 = title([' (t=',num2str(t(i)),')']);
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