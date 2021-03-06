% -----------------------------------------------------------------  %
% Matlab Programs included the Appendix B in the book:               %
%  Xin-She Yang, Engineering Optimization: An Introduction           %
%                with Metaheuristic Applications                     %
%  Published by John Wiley & Sons, USA, July 2010                    %
%  ISBN: 978-0-470-58246-6,   Hardcover, 347 pages                   %
% -----------------------------------------------------------------  %
% Citation detail:                                                   %
% X.-S. Yang, Engineering Optimization: An Introduction with         %
% Metaheuristic Application, Wiley, USA, (2010).                     %
%                                                                    % 
% http://www.wiley.com/WileyCDA/WileyTitle/productCd-0470582464.html % 
% http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0470582464.html  %
% -----------------------------------------------------------------  %
% ===== ftp://  ===== ftp://   ===== ftp:// =======================  %
% Matlab files ftp site at Wiley                                     %
% ftp://ftp.wiley.com/public/sci_tech_med/engineering_optimization   %
% ----------------------------------------------------------------   %


% Simulated Annealing for constrained optimization 
% by Xin-She Yang @ Cambridge University @2008
% Usage: sa_mincon(alpha)

function [bestsol,fmin,N,info]=sa_mincon_bldg(alpha,optParams,u_0)
ExactParams = [];
x0 = u_0;
%alpha = 0.8;
%make ub and lb for this fucker
Nu = size(optParams.B,2);
NL = optParams.len;

Lb = repmat(-optParams.U_feas.b(1:Nu),NL-1,Nu);
Ub = repmat(optParams.U_feas.b(Nu+1:end),NL-1,Nu);

Lb = Lb';
Ub = Ub';

u0 = x0;



% Default cooling factor
%if nargin<1, 
%    alpha=0.8; 
%end

% Display usage
disp('sa_mincon or [Best,fmin,N]=sa_mincon(0.8)');

% d dimensions

% % Welded beam design optimization
% Lb=[0.1 0.1  0.1 0.1];
% Ub=[2.0 10.0 10.0 2.0];
% u0=Lb+(Ub-Lb).*rand(size(Lb));

if length(Lb) ~=length(Ub),
    disp('Simple bounds/limits are improper!');
    return
end

%% Start of the main program -----------------------------------------
d=length(Lb);        % Dimension of the problem

% Initializing parameters and settings
eq_tol = 10^-2;
maxeval = 5000;
T_init = 1.0;        % Initial temperature
T_min =  1e-1;      % Finial stopping temperature
F_min = -1e+100;     % Min value of the function
max_rej=750;        % Maximum number of rejections
max_run=100;       % Maximum number of runs, 500 def
max_accept = 3000;   % Maximum number of accept, 250 def
initial_search=500; % Initial search period 
k = 5e-3;%1;               % Boltzmann constant
Enorm=1e-5;          % Energy norm (eg, Enorm=1e-8)

% variable for all accepts
CostSolutionList=Inf(maxeval,numel(u0)+1); %1st col for cost

% Initializing the counters i,j etc
i= 0; j = 0; accept = 0; totaleval = 0;
% Initializing various values
T = T_init;
E_init = Fun(u0,optParams,ExactParams);
E_old = E_init; E_new=E_old;
best=u0;  % initially guessed values
% Starting the simulated annealling
disp('...............Starting...............................');
tic
while (( (T > T_min) || (j <= max_rej)) && (totaleval<=maxeval))
    
    if(totaleval==0)
    disp('init cost')
    Fun(u0,optParams,ExactParams)
    
    end    
    i = i+1;
    
    if(rem(totaleval,100)==0)
       disp(num2str(totaleval)); 
    end
    % Check if max numbers of run/accept are met
    if (i >= max_run) || (accept >= max_accept)     
        % reset the counters
        i = 1;  accept = 1;
      % Cooling according to a cooling schedule
        T = cooling(alpha,T);  
        disp(strcat('The best found so far =',num2str(fmin)));
        if(fmin<eps)
           break; 
        end
    end
    
    % Function evaluations at new locations
    if totaleval<initial_search,
        init_flag=1;
        ns=newsolution(u0,Lb,Ub,init_flag);
    else
        init_flag=0;
        ns=newsolution(best,Lb,Ub,init_flag);
    end
     
      totaleval=totaleval+1;
      E_new = Fun(ns,optParams,ExactParams);
    % Decide to accept the new solution
    DeltaE=E_new-E_old;
    % Accept if improved
    if (DeltaE <0)
        best = ns; E_old = E_new;
        accept=accept+1;   j = 0;
        CostSolutionList(totaleval,:) = [E_new best];
       % disp('Accept w. improvement');
    else
%         j=j+1;
%         disp(num2str(DeltaE))
%         %keyboard
%         disp('Reject number')
%         j

    end
    % Accept with a small probability if not improved
    if (DeltaE>=0 && exp(-DeltaE/(k*T))>rand );
       best = ns; E_old = E_new;
       accept=accept+1;
       CostSolutionList(totaleval,:) = [E_new best];
       %disp('Accept w.o. improvement');
    else
       j=j+1;
       %disp(num2str(DeltaE))
       %keyboard
       %disp('Reject')
       %j
    end
    % Update the estimated optimal solution
    fmin=E_old;
    
end
info.CostSolutionList = CostSolutionList;
bestsol=best
fmin
N=totaleval
T
toc


%% New solutions
function s=newsolution(u0,Lb,Ub,init_flag)

  % Either search around
if length(Lb)>0 && init_flag==1,
  s=Lb+(Ub-Lb).*rand(size(u0));
else
  % Or local search by random walk
  s=u0+0.01*(Ub-Lb).*randn(size(u0));
end

s=bounds(s,Lb,Ub);

%% Cooling
function T=cooling(alpha,T)
T=alpha*T;

function ns=bounds(ns,Lb,Ub)
if length(Lb)>0,
% Apply the lower bound
  ns_tmp=ns;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
% Update this new move 
  ns=ns_tmp;
else
  ns=ns;
end


% d-dimensional objective function
function z=Fun(u,optParams,ExactParams)

% Objective
z=fobj(u,optParams,ExactParams);

% Apply nonlinear constraints by penalty method
% Z=f+sum_k=1^N lam_k g_k^2 *H(g_k) 
z=z+getnonlinear(u,optParams);

function Z=getnonlinear(u,optParams)
Z=0;
% Penalty constatn
lam=10^15; lameq=10^15;
[g,geq]=constraints(u,optParams);

% Inequality constraints
for k=1:length(g),
    Z=Z+ lam*g(k)^2*getH(g(k));
end

% Equality constraints (when geq=[], length->0)
for k=1:length(geq),
   Z=Z+lameq*geq(k)^2*geteqH(geq(k));
end


% Test if inequalities hold
function H=getH(g)
if g<=0, 
    H=0; 
else
    H=1; 
end

% Test if equalities hold
function H=geteqH(g)
eq_tol = 10^-2;
if abs(g)<=eq_tol,
    H=0;
else
    H=1; 
end


% Objective functions
function z=fobj(u,optParams,ExactParams)
% Welded Beam Design Optimization
% K. Ragsdell and D. Phillips, Optimal design of a class of welded
% strucures using geometric programming, 
% J. Eng. Ind., 98 (3):1021-1025, (1976).
% Best solution found in literature 
% [0.205730, 3.470489, 9.036624, 0.205729) with the objective 1.724852
% by Cagnina et al., Solving engineering optimization problems with the
% simple constrained particle swarm optimizer, Informatica, 32 (2008)319-326 

%z=1.10471*u(1)^2*u(2)+0.04811*u(3)*u(4)*(14.0+u(2));

%z = objfun_case_mex_onlne(u,optParams);
I1 = optParams.I1;
z = objfun_bldg(u,optParams,I1);
% All constraints
function [g,geq]=constraints(u,optParams)
[g,geq] = confun_bldg(u,optParams);
%% End of the program ------------------------------------------------

