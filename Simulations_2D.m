function Simulations_2D(in_K, in_model)
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
%%%  2D Simulations (Kelvin-Voigt and Maxwell models)                   %%%
%%%  For details about the equations and schemes please see the         %%%
%%%  manuscript indicated above and the supplementary material          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function Simulations_2D(in_K, in_model)
% Runs the spatially two-dimensional simulation using viscoelastic 
% model in_model (either 'Maxwell' or 'Kelvin-Voigt') on a uniform 
% square grid with in_K grid cells in each dimension. 
% The results are saved in file saved_y2D_[in_model]_[in_K].mat and 
% a video in saved_y2D_[in_model]_[in_K].avi .
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulations_2D.m: simulates a 2D mechanical model of pattern formation%
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
par.seed = 44;                   % Seed for random number generator
par.model = in_model;            % type of viscoelastic model
%par.model = 'Maxwell';
%par.model = 'Kelvin-Voigt'; 
par.K = in_K;                    % Number of spatial grid cells
par.L = 1;                       % Domain length
x1 = linspace(0,par.L,par.K+1);  % Discretise spatial domain
par.K2 = par.K^2;                % Total number of grid cells
par.dx = x1(2)-x1(1);            % Cell size
t0 = 0;                          % Initial time
tf = 10000;                      % Final time
tspan = linspace(t0,tf,251);     % Time span

%% Initial conditions
steadystate = [ones(2*par.K2,1); zeros(2*par.K2,1)];
rng(par.seed);                  % initialize random number generator                       % use seed 44 in random number generator
randP = [randn(par.K2,1); zeros(3*par.K2,1)]; % random perturbation
factP = 1e-2;                  % factor for random perturbation
y0 = steadystate+factP*randP;  % perturbed IC [n0, p0, u10, u20] as long column

%% Solve with ODE15i
solve_on = true;  % run simulation and save result: YES (true), NO (false)
if solve_on
    tic
    %%% Create Jacobian pattern
    [SPDY,SPDYP] = Jpattern(par.K);
    options = odeset('Jpattern',{SPDY,SPDYP});
    %options = [];
    %%% Compute consistent yp0
    %res = @(y,yp)(norm(mechanochemical(y,yp,par),'inf'));
    res = @(y,yp)(norm(mechanochemical(y,yp,par)));
    disp(['residuum of steady state = ' ...
        num2str(res(steadystate,0*steadystate), '%15.10e')]);
    [y0,yp0, resnorm] = decic(@(t,y,yp)(mechanochemical(y,yp,par)), ...
        t0, y0, ones(4*par.K2,1), ...
        zeros(4*par.K2,1), zeros(4*par.K2,1), options);
    disp(['residuum (from decic) of IC = ' num2str(resnorm, '%15.10e')]);
    disp(['residuum (from res()) of IC = ' num2str(res(y0,yp0), '%15.10e')]);
    toc
    %%% Solve
    tic
    [t,y] = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,...
        y0,yp0,options);
    toc
    %%% Save computed solution to file
    filename = ['saved_y2D_' par.model '_' num2str(par.K)];
    save(filename, 't', 'y', 'par', 'x1');
end


%%% Plot
if isempty(filename)
  filename = 'saved_y2D_Maxwell_50'; % make sure it exists...
end
load(filename, 't', 'y', 'par', 'x1');
video_on = true; % Record video: YES (true), NO (false)
video_filename = [filename '.avi'];
plot_solution(x1,y,t,par,video_on,video_filename);

end

%% Main function implementing the model
function f = mechanochemical(y,yp,par)

    %%% Parameter values
    eta = 1;           % viscosity
    E = 1;             % elasticity
    D = 0.01;          % diffusion
    alpha = 0.05;      % haptotaxis
    r = 1;             % proliferation
    s = 10;            % substrate elasticity
    beta = 0.005;      % long range traction 
    lambda = 0.5;      % (cell traction) saturation coefficient
    tau = 0.2;         % cell traction 
    nu = 0.25;         % Poisson ratio
    Ep = E/(1+nu);     % eq.(33)
    nup = nu/(1-2*nu); % eq.(33)
    % mu = nup*eta:    % Simplifying assumption (A.4)

    %%% Choose constitutive model - see eq.(31) and Table 3
    switch par.model
        case 'Kelvin-Voigt' % Kelvin Voigt  - (A.2)
            [a0,a1,b0,b1,c0,c1] = deal(1/eta,0,Ep/eta,1,Ep*nup/eta,nup); 
        case 'Maxwell' % Maxwell - (A.3)
            [a0,a1,b0,b1,c0,c1] = deal(1/eta,1/Ep,0,1,0,nup);
        otherwise
            error('Unknown constitutive model')
    end
    %%% Reshape input vectors
    [n,p,u1,u2] = vec2mats(y,par.K);
    [np,pp,u1p,u2p] = vec2mats(yp,par.K);
    
    %%% Equation for n
    % Advection velocity at grid cell interfaces - eq.(S.22)
    u1p_Avx1 = Avx1(u1p);
    vx1 = alpha*Mx1_face(p,par)+u1p_Avx1; 
    u2p_Avx2 = Avx1(u2p')';
    vx2 = alpha*Mx1_face(p',par)'+u2p_Avx2;
    % fn(t,n,n') = 0 - eq.(S.21)
    fn = np ...
        - D*(Mxx1(n,par) + Mxx1(n',par)') ...
        +  MA1(n, vx1, par) + MA1(n', vx2', par)' ...
        - r * (n.*(1-n));
      
    %%% Equation for p
    % fp(t,p,p') = 0 - eq.(S.23)
    fp = pp ...
        + MA1(p, u1p_Avx1, par) + MA1(p', u2p_Avx2', par)';
    
    %%% Equation for u1 and u2    
    % Traction term - eq.(S.29) with (S.13) and (S.30)
    pexp = 2;
    fn1 = n./(1+lambda*n.^pexp); % Lambda_1
    fn2 = ((1-(pexp-1)*lambda*n.^pexp) ./ ((1+lambda*n.^pexp).^2)); % Lambda_2 
    fp1 = p + beta*(Mxx1(p, par)+ Mxx1(p', par)'); % M_T2 P
    fp2 = pp + beta*(Mxx1(pp, par)+ Mxx1(pp', par)');  % M_T2 P'
    Tr = tau*(a0*fn1.*fp1 + a1*(fn2.*np.*fp1+fn1.*fp2));    
    % fu1(n,n',p,p',u1,u1',u2,u2') = 0 - eq.(S.27)
    fu1 = (b1+c1)*Mxx1(u1p, par) + 0.5*b1*(Mxx1(u1p', par))' ...
        +(0.5*b1+c1)* Mx1x2(u2p, par) ...
        +(b0+c0)*Mxx1(u1, par)   + 0.5*b0*(Mxx1(u1', par))' ...
        +(0.5*b0+c0)* Mx1x2(u2, par) ...
        + Mx1_center(Tr, par) ...
        -a1*s*(p.*u1p + pp.*u1) - a0*s*p.*u1;
    % fu2(n,n',p,p',u1,u1',u2,u2') = 0 - eq.(S.28)
    fu2 = (b1+c1)*Mxx1(u2p', par)' + 0.5*b1*Mxx1(u2p, par) ...
        +(0.5*b1+c1)* Mx1x2(u1p, par) ...
        +(b0+c0)*Mxx1(u2', par)'   + 0.5*b0*Mxx1(u2, par) ...
        +(0.5*b0+c0)* Mx1x2(u1, par) ...
        +Mx1_center(Tr', par)'  ...
        -a1*s*(p.*u2p + pp.*u2) - a0*s*p.*u2;
    
    %%% Full system - eq.(S.15)
    f = mats2vec(fn, fp, fu1, fu2);
    
end
    
%% Annexed functions

%%% Vector to matrices
function [n,p,u1,u2] = vec2mats(y,K)
    K2 = K^2;
    n  = reshape(y(0*K2+1:1*K2), K, K);
    p  = reshape(y(1*K2+1:2*K2), K, K);
    u1 = reshape(y(2*K2+1:3*K2), K, K);
    u2 = reshape(y(3*K2+1:4*K2), K, K);
end

%%% Matrices to vector
function y = mats2vec(n, p, u1, u2)
    y = [n(:); p(:); u1(:); u2(:)];
end

%%% Compute advection at grid cell centers using first order upwinding
%%% in the x1 direction with advective velocity given at grid cell 
%%% interfaces - def.(S.7) and (S.20)
function fluxdiffx1 = MA1(y, vel, par)
    % find pos/neg velocity entries
    velPos_ind = (vel>0);
    velNeg_ind = (~velPos_ind);
    % compute flux accross cell interfaces using first order upwinding 
    % def.(S.8)-(S.9)
    flux = NaN(size(y));
    flux(velPos_ind) = y(velPos_ind).*vel(velPos_ind);
    yshift = y([2:end 1],:);
    flux(velNeg_ind) = yshift(velNeg_ind).*vel(velNeg_ind);
    % add entries for flux along x1=0 by periodicity
    flux = flux([end,1:end],:);
    % compute flux difference per grid cell in x1 direction - def.(S.7) and (S.4)
    fluxdiffx1 = (1/par.dx)*(flux(2:end,:)-flux(1:(end-1),:));
end

%%% Compute second order mixed derivative in grid cell centres
%%% (second order approximation, central) - def.(S.18)
function dx1x2 = Mx1x2(y, par)
    persistent Mdx;
    if size(Mdx,2) ~= size(y,1)
      disp('Mx1x2     : resetting persistent matrix.');
      c = (1/(2*par.dx))*sparse([0;-1; zeros(par.K-3,1); 1]);
      Mdx = toeplitz(c,[c(1), c(end:-1:2)']); 
    end
    dx1x2 = (Mdx*(Mdx*y)')';
end

%%% Compute first order derivative in x1 direction on grid cell interfaces
%%% (second order approximation, central) - def.(S.16)
function dx1 = Mx1_face(y, par)
    persistent Mdx;
    if size(Mdx,2) ~= size(y,1)
      disp('Mx1_face  : resetting persistent matrix.');
      c = (1/par.dx)*sparse([-1; zeros(par.K-2,1); 1]);
      Mdx = toeplitz(c,[c(1), c(end:-1:2)']); 
    end
    dx1 = Mdx*y;
end

%%% Compute first order derivative in x1 direction on grid cell centres
%%% (second order approximation, central) - def.(S.16)
function dx1 = Mx1_center(y, par)
    persistent Mdx;
    if size(Mdx,2) ~= size(y,1)
      disp('Mx1_center: resetting persistent matrix.');
      c = (1/(2*par.dx))*sparse([0;-1; zeros(par.K-3,1); 1]);
      Mdx = toeplitz(c,[c(1), c(end:-1:2)']); 
    end
    dx1 = Mdx*y;
end

%%% Compute variable at grid cell interfaces in the x1 direction -
%%% def.(S.3) and (S.19)
function avx1 = Avx1(y)
    avx1 = 0.5*(y +y([2:end,1],:));
end

%%% Compute second order derivative in x1 direction in grid cell centres
%%% (second order approximation, central) - def.(S.2) and (S.17)
function dxx1 = Mxx1(y, par)
    persistent Mdxx;
    if size(Mdxx,2) ~= size(y,1)
      disp('Mxx1      : resetting persistent matrix.');
      c = (1/(par.dx^2))*sparse([-2; 1; zeros(par.K-3,1); 1]);
      Mdxx = toeplitz(c,[c(1), c(end:-1:2)']);
    end
    dxx1 = Mdxx*y;    
end

%%% Plot solution 
function plot_solution(x1,y,t,par,video_on,video_filename)
    if video_on % Initialise video
        vid = VideoWriter(video_filename);
        open(vid);
        figure('Units','normalized','Position',[0 0 0.5 0.45])
    end
    for i=1:length(t)
        clf
        [N, P, U1, U2] = vec2mats(y(i,:),par.K);
        n = [N(par.K,par.K)  N(par.K,:)
             N(:,par.K)      N];
        p = [P(par.K,par.K)  P(par.K,:)
             P(:,par.K)      P];
        u1 = [U1(par.K,par.K) U1(par.K,:)
              U1(:,par.K)     U1];
        u2 = [U2(par.K,par.K) U2(par.K,:)
              U2(:,par.K)     U2];
        subplot(2,2,1)
        surf(x1,x1,n)
        view(0,90)
        shading flat
        axis square
        xlim([0,par.L])
        ylim([0,par.L])
        title('$n(t,x)$')
        colormap(parula);
        colorbar;
        subplot(2,2,2)
        surf(x1,x1,p)
        view(0,90)
        shading flat
        axis square
        xlim([0,par.L])
        ylim([0,par.L])
        title('$\rho(t,x)$')
        colormap(parula);
        colorbar;
        subplot(2,2,3)
        surf(x1,x1,u1)
        view(0,90)
        shading flat
        axis square
        xlim([0,par.L])
        ylim([0,par.L])
        title('$u_1(t,x)$')
        colormap(parula);
        colorbar;
        subplot(2,2,4)
        surf(x1,x1,u2)
        view(0,90)
        shading flat
        axis square
        xlim([0,par.L])
        ylim([0,par.L])
        title('$u_2(t,x)$')
        colormap(parula);
        colorbar;
        a = axes;
        t1 = title([par.model ' (t=',num2str(t(i)),')'],...
            'Position', [0.5, 1, 2]);
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

function [SPDY,SPDYP] = Jpattern(K)
    % init return matrices to empty sparse matrices
    SPDY  = sparse(4*K^2,4*K^2); 
    SPDYP = sparse(4*K^2,4*K^2);
    
    % vector of first entry (start index) of n, p, u1, u2 in y and yp
    si = 1+(0:3)*(K^2);
    % vector of the "grid cell numbers-1" from 0 to K^2-1; this vector
    % represents the order of the grid cells
    self = (0:(K^2-1))'; 
    % indices of all entries in y and yp of the four solution components 
    in=si(1)+self;
    ip=si(2)+self;
    iu1=si(3)+self;
    iu2=si(4)+self;
   
%     % Example stencil: second order central difference of 2D Laplace
%     centre = [2,2]; % [row index, column index] of stencil centre
%     % stencil: visual mode, i.e. increasing row index is decreasing
%     % x2-direction and increasing column index is increasing x1-direction
%     stencil = [0 1 0
%                1 1 1
%                0 1 0];
           
    % Here: for now generic full stencil of size 5x5
    stencilsize = 5;
    centre = [1, 1]*((stencilsize+1)/2);
    stencil = ones(stencilsize); % square matrix of all 1
    
    % rather simple, not optimal but sufficient choice: 
    % assume each equation depends on n, p, u1, u2 and their 
    % derivatives using the same stencil as defined above.
    mat = SPDY;
    mat = add(mat, in, in, stencil, centre, K);
    mat = add(mat, in, ip, stencil, centre, K);
    mat = add(mat, in, iu1, stencil, centre, K);
    mat = add(mat, in, iu2, stencil, centre, K);
    %
    mat = add(mat, ip, in, stencil, centre, K);
    mat = add(mat, ip, ip, stencil, centre, K);
    mat = add(mat, ip, iu1, stencil, centre, K);
    mat = add(mat, ip, iu2, stencil, centre, K);
    %
    mat = add(mat, iu1, in, stencil, centre, K);
    mat = add(mat, iu1, ip, stencil, centre, K);
    mat = add(mat, iu1, iu1, stencil, centre, K);
    mat = add(mat, iu1, iu2, stencil, centre, K);
    %
    mat = add(mat, iu2, in, stencil, centre, K);
    mat = add(mat, iu2, ip, stencil, centre, K);
    mat = add(mat, iu2, iu1, stencil, centre, K);
    mat = add(mat, iu2, iu2, stencil, centre, K);
    SPDY=mat;
    
    mat = SPDYP;
    mat = add(mat, in, in, stencil, centre, K);
    mat = add(mat, in, ip, stencil, centre, K);
    mat = add(mat, in, iu1, stencil, centre, K);
    mat = add(mat, in, iu2, stencil, centre, K);
    %
    mat = add(mat, ip, in, stencil, centre, K);
    mat = add(mat, ip, ip, stencil, centre, K);
    mat = add(mat, ip, iu1, stencil, centre, K);
    mat = add(mat, ip, iu2, stencil, centre, K);
    %
    mat = add(mat, iu1, in, stencil, centre, K);
    mat = add(mat, iu1, ip, stencil, centre, K);
    mat = add(mat, iu1, iu1, stencil, centre, K);
    mat = add(mat, iu1, iu2, stencil, centre, K);
    %
    mat = add(mat, iu2, in, stencil, centre, K);
    mat = add(mat, iu2, ip, stencil, centre, K);
    mat = add(mat, iu2, iu1, stencil, centre, K);
    mat = add(mat, iu2, iu2, stencil, centre, K);
    SPDYP=mat;
    
    % final tuning
    SPDYP(SPDYP>1) = 1;
    SPDY(SPDY>1) = 1;
   
    function M = add(M, c, r, stencil, centre, K)
        toR = @(ind,count)(reshape(circshift(reshape(ind,K,K),-count,1),1,K^2));
        toU = @(ind,count)(reshape(circshift(reshape(ind,K,K),-count,2),1,K^2));
        [x2dir,x1dir]= find(stencil);
        x1dir = x1dir-centre(2);
        x2dir = -(x2dir-centre(1));
        %[x1dir,x2dir]
        for ii=1:length(x1dir)
            rtmp = toR(r   , x1dir(ii));
            rtmp = toU(rtmp, x2dir(ii));
            M = M + sparse(c,rtmp,ones(size(c)),4*K^2,4*K^2);
        end
    end
end
