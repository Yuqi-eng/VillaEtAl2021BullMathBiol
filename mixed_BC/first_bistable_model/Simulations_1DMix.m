function Simulations_1DMix(in_K)
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
tf = 1000;                        % Final time
tspan = linspace(t0,tf,100);     % Time span
dt = tspan(2)-tspan(1);

%% Initial conditions - eq.(28)
steadystate = [ones(2*par.K,1); zeros(par.K,1)];
f = 1+5*exp(-((par.x)./0.2).^2);
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
% res = @(y,yp)(norm(mechanochemical(y,yp,par)));
% disp(['residuum of steady state = ' ...
%     num2str(res(steadystate,0*steadystate), '%15.10e')]);
% opt = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7));
% nfixed = zeros(par.K,1);
% pfixed = zeros(par.K,1);
% ufixed = zeros(par.K,1);
% [y0,yp0,resnorm] = decic(@(t,y,yp)(mechanochemical(y,yp,par)), t0, ...
%     y0, [nfixed; pfixed; ufixed], ...
%     [zeros(3*par.K,1)], zeros(3*par.K,1), opt);
% disp(['residuum (from decic) of IC = ' num2str(resnorm, '%15.10e')]);
% disp(['residuum (from res()) of IC = ' num2str(res(y0,yp0), '%15.10e')]);

%% Solve
tic
% [t,y] = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
a = [1 2 3 4 5];
tau1 = [0.1 0.2 0.3 0.4 0.5];
nleft = [];
for i = 1:5
    for j = 1:5
    par.a = a(i);
    par.tau1 = tau1(j);
    %% Solve with ODE15i
    %%% Compute consistent yp0 
    res = @(y,yp)(norm(mechanochemical(y,yp,par)));
    disp(['residuum of steady state = ' ...
        num2str(res(steadystate,0*steadystate), '%15.10e')]);
    opt = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7));
    nfixed = zeros(par.K,1);
    pfixed = zeros(par.K,1);
    ufixed = zeros(par.K,1);
    [y0,yp0,resnorm] = decic(@(t,y,yp)(mechanochemical(y,yp,par)), t0, ...
        y0, [nfixed; pfixed; ufixed], ...
        [zeros(3*par.K,1)], zeros(3*par.K,1), opt);
    disp(['residuum (from decic) of IC = ' num2str(resnorm, '%15.10e')]);
    disp(['residuum (from res()) of IC = ' num2str(res(y0,yp0), '%15.10e')]);
    sol = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
    tfinal = sol.stats.tfinal;
    t = linspace(t0,tfinal,tfinal/dt+1);
    [y,yp] = deval(sol,t);
    nleft(end+1) = y(1,end);
    end
end
disp(nleft);
toc
%%% Save computed solution to file
filename = ['saved_y1D_' num2str(par.K)];
save(filename, 't', 'y', 'par', 'x', 'yp');

%% Plot
% video_on = true; % Record video: YES (true), NO (false)
% video_filename = [filename '.avi'];
% plot_solution(x,y,yp,t,par,video_on,video_filename);
pic_name = ['res' '.png'];
% plot_res(x,y,yp,par,tf,pic_name);

end

%% Main function implementing the model
function f = mechanochemical(y,yp,par)
    %%% Parameter values
    a = par.a;
    D2 = 0.1;
    d2 = 1;
    E1 = 1.4e8;
    tau1 = par.tau1;
    
    %%% Reshape input vectors
    [n,p,u] = deal(y(1:par.K),y(par.K+1:2*par.K),y(2*par.K+1:3*par.K));
    ntilde = [n(2); n; par.Nright]; 
    ptilde = [p(2); p; par.Pright]; 
    utilde = [0; u; 0];
    [np,pp,up] = deal(yp(1:par.K),yp(par.K+1:2*par.K),...
         yp(2*par.K+1:3*par.K));
    uptilde = [0; up; 0];

    %%% Traction force term
    % Tr = tau*p.*n;
    % Trtilde = [Tr(2); Tr; tau*ptilde(end)*ntilde(end)];
    % for when traction force term has a hill function
    n0 = 1.25;
    k2 = 4;
    h1 = (n.^k2)./(n0^k2*ones(size(n)) + n.^k2);
    Tr = tau1.*p.*h1;
    n2 = (ntilde(end)^k2)/(n0^k2 + ntilde(end)^k2);
    Trtilde = [Tr(2); Tr; tau1*ptilde(end)*n2];

    %%% Equation for n
    sig = Mx(uptilde, par) + E1*Mx(utilde,par) + Tr;
    % disp(sig(1));
    sig0 = 0.3*ones(size(sig));
    k1 = 4;
    fsig = (sig.^k1)./(sig0.^k1 + sig.^k1);
    % disp(fsig.');
    fn = np - Mxx(ntilde, par) + MA(n,up,par) - a*fsig - ones(size(n)) + n;

    %%% Equation for p
    fp = pp - D2*Mxx(ptilde,par) + MA(p,up,par) - n + d2*p;

    %%% Equation for u 
    fu = Mxx(uptilde, par) + E1*Mxx(utilde,par) + tau1*Mx(Trtilde, par);

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

%%% Compute flux gradient using first order upwinding without center
%%% averaging
%%% Take material at K grid cells
function fluxdiffx1 = MA(y, vel, par)
    flux = zeros(size(y));

    for i=2:size(flux)
        if vel(i)>0
            flux(i) = y(i-1)*vel(i);
        else
            flux(i) = y(i)*vel(i);
        end
    end

    % flux(1) = 0;
    flux(1) = vel(1)*y(1);

    % compute flux difference per grid cell - def.(S.7) and (S.4)
    fluxdiffx1 = Mx([0; flux; 0],par);
    fluxdiffx1(1) = (flux(2)-flux(1))/par.dx;
    fluxdiffx1(end) = (flux(end)-flux(end-1))/par.dx;

    % vx = Mx([0;vel;0],par);
    % fluxdiffx1(end) = vx(end)*1;

    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is 0
    % fluxdiffx1 = [(1/par.dx)*(flux(2:end)-flux(1:(end-1))); 0];

    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is the same as to its left 
    % fluxdiffx1 = (1/par.dx)*(flux(2:end)-flux(1:(end-1)));
    % fluxdiffx1 = [fluxdiffx1; fluxdiffx1(end)];
end

%%% Plot solution 
function plot_solution(x,y,yp,t,par,video_on,video_filename)
    if video_on % Initialise video
        vid = VideoWriter(video_filename);
        open(vid);
        figure('Units','normalized','Position',[0 0 0.5 0.45])
    end
    maxu = 10^(-6);

    for i=1:length(t)
        clf
        n = [y(1:par.K,i)];
        p = [y(par.K+1:2*par.K,i)];
        u = [y(2*par.K+1:3*par.K,i)];
        v = [yp(2*par.K+1:3*par.K,i)];
        if max(abs(u))>maxu
          maxu = max(abs(u))+0.1*max(abs(u));
        end
        subplot(1,4,1)
        plot(x,n)
        title('$n(t,x)$')
        % ylim([0.8,max(1.2,max(n))])
        axis square
        subplot(1,4,2)
        plot(x,p)
        title('$\rho(t,x)$')
        % ylim([0.8,max(1.2,max(p))])
        axis square
        subplot(1,4,3)
        plot(x,u)
        title('$u(t,x)$')
        axis square
        subplot(1,4,4)
        plot(x,v)
        title('$v(t,x)$')
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

%%% Plot solution 
function plot_res(x,y,yp,par,tf,pic_name)
    figure('Units','normalized','Position',[0 0 0.5 0.45])
    n = [y(1:par.K,end)];
    p = [y(par.K+1:2*par.K,end)];
    u = [y(2*par.K+1:3*par.K,end)];
    v = [yp(2*par.K+1:3*par.K,end)];
    % res_file = 'res';
    % save(res_file,"n","p","u");
    subplot(1,3,1)
    plot(x,n)
    % title('$n(t,x)$')
    title('Cell Density')
    ylim([0,3])
    axis square
    subplot(1,3,2)
    plot(x,p)
    % title('$\rho(t,x)$')
    title('Collagen Density')
    ylim([0,3])
    axis square
    subplot(1,3,3)
    plot(x,u)
    % title('$u(t,x)$')
    title('Displacement field')
    ylim([-1,1])
    axis square
    % subplot(1,4,4)
    % plot(x,v)
    % title('$v(t,x)$')
    % axis square
    a = axes;
    a.Visible = 'off';
    t1 = title([' (t=',num2str(tf),')']);
    t1.Visible = "on";
    saveas(gcf,pic_name);
end