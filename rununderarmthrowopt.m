
%% Underarm Throwing simulation
% Performs analysis of the underarm throwing model by Alexander,
% including different types of optimization. The first uses Alexander's
% original proximal-to-distal activation sequence, in which the shoulder
% is activated maximally, and the elbow is also activated maximally after
% a delay relative to the shoulder. The second uses a piecewise linear
% spline for the elbow, allowing for arbitrary activation values over
% time; the shoulder is kept maximal.

function rununderarmthrowopt(whichparts)
% set whichparts = [1 2 3 ...] to set which parts to run (default all):
%   1 Alexander's proximal-to-distal sequence, maximal 
%   2 Piecewise linear spline for elbow
%   3 Cubic spline for elbow
%   4 Piecewise linear splines for shoulder and elbow
%   5 Cubic splines for shoulder and elbow

bestthrows = NaN*ones(1,5); % store max throw velocities here
timings = NaN*ones(1,5);    % store timing info here
if nargin == 0
  whichparts = 1:5;
end
%% Model parameters
% These parameters determine the geometry and inertia of the arm
% segments, and are stored in a structure to be passed to the
% simulation code. Also included here are some simulation parameters,
% for max time to simulate and the initial state vector.

% arm segment length, mass
s = 1; M = 1; g = 1; % we are neglecting gravity, but g is included for 
% non-dimensionalization

% store these in the parameters stucture
armparms.s = s; armparms.M = M; armparms.g = g; 
% or, to use real units, just set these values:
% armparms.s = 0.44; armparms.M = 70; armparms.g = 9.81;

armparms.G = 3; % G is force-velocity shape parameter

% thkdmax is the max angular velocity of knee, thedmax of elbow
% describing force-velocity curve in angular velocities
armparms.thsdmax = 8*sqrt(g/s); armparms.thedmax = 8*sqrt(g/s);

% Max torque at knee and elbow, given in terms of M*g*s
armparms.Tsmax = 0.7*M*g*s; armparms.Temax = 0.7*M*g*s;

armparms.lua = 0.8*s; armparms.lfa = 0.8*s;
armparms.rua = 0.41*s; armparms.rfa = 0.28*s;
%rua = 0.4*s; rfa = 0.8*s;
armparms.mua = 0.033*M; armparms.mfa = 0.021*M; armparms.mball = 0.011*M;

% The following are parameters related to the simulation itself.
armparms.maxtime = 1.4;   % maximum simulation time
armparms.x0 = [0 0 0 0];  % initial state for the simulation
% where state is [qshoulder qelbow qdotshoulder qdotelbow]
 
%% Alexander's proximal-to-distal activation sequence
% Uses maximally activated shoulder and
% elbow, with elbow activation delayed by time tdelay. These are
% simulated as two epochs, first using a 1-dof model with the two
% segments sharing the same acceleration, effectively keeping the
% elbow straight. Afer tdelay the simulation is stopped and the final
% conditions used as initial conditions to the second epoch, where
% both torques are activated and the segments simulated as a 2-dof
% system.

% Compute a parameter sweep of simulations with different delays
% and find the maximum ball speed from each
if any(whichparts == 1) 
tic
Tdelays = linspace(0.02, 1, 30); % Beyond delay of 1, shoulder goes too far
for i = 1:length(Tdelays)
  Tdelay = Tdelays(i);
  maxvs(i) = -underarmthrowalexander(Tdelay, armparms, 0);
end

% Find the optimum throw using fminsearch:
%
% Define an anonymous function to allow fminsearch to find an optimum
optfunction = @(optarg) underarmthrowalexander(optarg, armparms, 0);
% (thus optfunction(optvec) is equivalent to
% underarmthrowalexander(optvec, armparms, 0)

optvec = 0.8; % initial guess a delay of 0.8
[optstar, vstar, exitflag] = ...
  fminsearch(optfunction, optvec);

figure(1); clf;
plot(Tdelays, maxvs, '-', optstar, -vstar, '*');
legend('Max velocity', ...
  sprintf('vstar = %4.2f  delay = %4.2f', -vstar, optstar));
xlabel('delay time'); ylabel('max ball velocity')
snapnow

% Simulate the optimum and plot the results
figure(2); clf;
[negmaxv, alex.ts, alex.xs] = underarmthrowalexander(optstar, armparms, 1);
alex.maxv = -negmaxv; snapnow;

% Also show an animation on the screen (but not in the publish document)
figure(3); clf;
% the simulation includes many time points, which need to be down-
% sampled for the animation or cartoon
[tsinterp, iinterp] = unique(alex.ts); xsinterp = alex.xs(iinterp,:);
xsdownsampled = interp1(tsinterp,xsinterp,linspace(0,tsinterp(end),14));
animatemodel(xsdownsampled, armparms);

% And draw a stop motion cartoon of the optimum
clf; % where we re-sample the states to give a reasonable cartoon
cartoonmodel(xsdownsampled, armparms);
title(sprintf('Maximal throw vopt = %4.2f', alex.maxv));
bestthrows(1) = alex.maxv;
timings(1) = toc;
end % alexander part

%% Simulation with piecewise linear spline for elbow
% Keep the shoulder maximally activated throughout, but for the elbow,
% define a series of knot points, which are interpolated between with
% piecewise linear curves. This allows for more flexibility than the
% delayed elbow activation.

if any(whichparts ==2)
tic

% Choose how many knots to include
nknots = 10; linearsplineparms.nknots = nknots;
linearsplineparms.splinetype = 2; % spline order is this number minus 1
linearsplineparms.timevec = linspace(0, armparms.maxtime, linearsplineparms.nknots);

% Set up the initial guess at some knot points
shoulderknots = armparms.Tsmax * ones(1,nknots); % all maximal
elbowknots = armparms.Temax * ones(1,nknots);    % all maximal, except
% start with low force for the first part of the throw
elbowknots(linearsplineparms.timevec < 0.8) = 0.15*armparms.Temax; 

% Perform an optimization on the elbowknots to find the max velocity.
% The code below sets up fminsearch.

optvec = [elbowknots]; % this is the vector of values to be optimized

% to set up the optimization with fmincon, define lower and upper
% bounds on the optvec values:
lowerbound = optvec * 0; upperbound = ones(1,nknots)*armparms.Temax;

options = optimoptions(@fmincon, 'Display', 'iter', 'MaxFunEvals', 2500, ...
  'Algorithm', 'active-set');

% The optimizing function is underarmthrowsplineelbow, called through an
% anonymous function so that the parameters are passed to it.
optfunction = @(optarg) ...
  underarmthrowsplineelbow(optarg, armparms, linearsplineparms, 0);
% (thus optfunction(optvec) is equivalent to
% underarmthrowsplineelbow(optvec, armparms, linearsplineparms, 0)

% Here is the optimization using fmincon (constrained optimization)
[optstar, vstar, exitflag, output] = fmincon(optfunction, optvec, [], ...
  [], [], [], lowerbound, upperbound, [], options);

figure(1); clf;
[negmaxv, elbow.ts, elbow.xs] = ...
  underarmthrowsplineelbow(optstar, armparms, linearsplineparms, 1);
% note that the results ts and xs contain time and states only for the
% knot points, so there's no need to down-sample for animation
elbow.maxv = -negmaxv;
% draw the animation or stop motion cartoon
figure(2); clf;
%animatemodel(xs, armparms); % uncomment if you want to see animation
cartoonmodel(elbow.xs, armparms);
title(sprintf('Maximal throw vopt = %4.2f', elbow.maxv));
bestthrows(2) = elbow.maxv;
timings(2) = toc;
end % elbow piecewise linear part

%% Simulation with cubic spline for elbow
% Keep the shoulder maximally activated throughout, but for the elbow,
% define a series of knot points, which are interpolated between with
% cubic curves. This allows for more flexibility than the
% delayed elbow activation.

if any(whichparts == 3)
tic
% Choose how many knots to include
nknots = 10; cubicsplineparms.nknots = nknots;
cubicsplineparms.splinetype = 4; % cubic spline requires value of 4
cubicsplineparms.timevec = linspace(0, armparms.maxtime, cubicsplineparms.nknots);

% Set up the initial guess at some knot points
shoulderknots = armparms.Tsmax * ones(1,nknots); % all maximal
%elbowknots = optstar; % use knots from piecewise linear spline as guess
elbowknots = armparms.Temax * ones(1,nknots);    % all maximal, except
% start with low force for the first part of the throw
elbowknots(cubicsplineparms.timevec < 0.8) = 0.15*armparms.Temax; 

% Perform an optimization on the elbowknots to find the max velocity.
% The code below sets up fminsearch.

optvec = [elbowknots]; % this is the vector of values to be optimized

% to set up the optimization with fmincon, define lower and upper
% bounds on the optvec values:
lowerbound = optvec * 0; upperbound = ones(1,nknots)*armparms.Temax;

options = optimoptions(@fmincon, 'Display', 'iter', 'MaxFunEvals', 2500, ...
  'Algorithm', 'active-set');

% The optimizing function is underarmthrowsplineelbow, called through an
% anonymous function so that the parameters are passed to it.
optfunction = @(optarg) ...
  underarmthrowsplineelbow(optarg, armparms, cubicsplineparms, 0);
% (thus optfunction(optvec) is equivalent to
% underarmthrowsplineelbow(optvec, armparms, cubicsplineparms, 0)

% Here is the optimization using fmincon (constrained optimization)
[optstar, vstar, exitflag, output] = fmincon(optfunction, optvec, [], ...
  [], [], [], lowerbound, upperbound, [], options);

figure(1); clf;
[negmaxv, elbow.ts, elbow.xs] = ...
  underarmthrowsplineelbow(optstar, armparms, cubicsplineparms, 1);
% note that the results ts and xs contain time and states only for the
% knot points, so there's no need to down-sample for animation
cubicelbow.maxv = -negmaxv;
% draw the animation or stop motion cartoon
figure(2); clf;
%animatemodel(xs, armparms); % uncomment if you want to see animation
cartoonmodel(elbow.xs, armparms);
title(sprintf('Maximal throw vopt = %4.2f', cubicelbow.maxv));
bestthrows(3) = cubicelbow.maxv;
timings(3) = toc;

end % cubic spline elbow part

%% Simulation with double piecewise linear splines
% Both shoulder and elbow controlled by piecewise linear trajectories.
% Define a series of knot points, which are interpolated between with
% piecewise linear curves. This allows for more flexibility than the
% singe elbow spline.

if any(whichparts == 4)
tic
% Choose how many knots to include
nknots = 10; linearsplineparms.nknotstot = 2*nknots; % two sets of knots
linearsplineparms.splinetype = 2; % spline order is this number minus 1
linearsplineparms.timevec = linspace(0, armparms.maxtime, nknots);

% Set up the initial guess at some knot points
shoulderknots = armparms.Tsmax * ones(1,nknots); % all maximal
elbowknots = armparms.Temax * ones(1,nknots);    % all maximal, except
% start with low force for the first part of the throw
elbowknots(linearsplineparms.timevec < 0.8) = 0.15*armparms.Temax; 

% Perform an optimization on the elbowknots to find the max velocity.
% The code below sets up fminsearch.

optvec = [shoulderknots elbowknots]; % this is the vector of values to be optimized

% to set up the optimization with fmincon, define lower and upper
% bounds on the optvec values:
lowerbound = optvec * 0; 
upperbound = [ones(1,nknots)*armparms.Tsmax ones(1,nknots)*armparms.Temax];

options = optimoptions(@fmincon, 'Display', 'iter', 'MaxFunEvals', 2500, ...
  'Algorithm', 'active-set');

% The optimizing function is underarmthrowsplineboth, called through an
% anonymous function so that the parameters are passed to it.
optfunction = @(optarg) ...
  underarmthrowsplineboth(optarg, armparms, linearsplineparms, 0);

% Here is the optimization using fmincon (constrained optimization)
[optstar, vstar, exitflag, output] = fmincon(optfunction, optvec, [], ...
  [], [], [], lowerbound, upperbound, [], options);

figure(1); clf;
[negmaxv, both.ts, both.xs] = ...
  underarmthrowsplineboth(optstar, armparms, linearsplineparms, 1);
% note that the results ts and xs contain time and states only for the
% knot points, so there's no need to down-sample for animation
both.maxv = -negmaxv;
% draw the animation or stop motion cartoon
figure(2); clf;
%animatemodel(xs, armparms); % uncomment if you want to see animation
cartoonmodel(both.xs, armparms);
title(sprintf('Maximal throw vopt = %4.2f', both.maxv));
bestthrows(4) = both.maxv;
timings(4) = toc;
end % double piecewise linear

%% Simulation with double cubic splines
% Both shoulder and elbow controlled by cubic spline trajectories.
% Define a series of knot points, which are interpolated between with
% cubic spline curves. This allows for more flexibility than the
% singe elbow spline.

if any(whichparts == 5)
tic
% Choose how many knots to include
nknots = 10; cubicsplineparms.nknotstot = 2*nknots; % two sets of knots
cubicsplineparms.splinetype = 4; % spline order is this number minus 1
cubicsplineparms.timevec = linspace(0, armparms.maxtime, nknots);

% Set up the initial guess at some knot points
shoulderknots = armparms.Tsmax * ones(1,nknots); % all maximal
elbowknots = armparms.Temax * ones(1,nknots);    % all maximal, except
% start with low force for the first part of the throw
elbowknots(cubicsplineparms.timevec < 0.8) = 0.15*armparms.Temax; 

% Perform an optimization on the spline knots to find the max velocity.
optvec = [shoulderknots elbowknots]; % this is the vector of values to be optimized

% to set up the optimization with fmincon, define lower and upper
% bounds on the optvec values:
lowerbound = optvec * 0; 
upperbound = [ones(1,nknots)*armparms.Tsmax ones(1,nknots)*armparms.Temax];

options = optimoptions(@fmincon, 'Display', 'iter', 'MaxFunEvals', 2500, ...
  'Algorithm', 'active-set');

% The optimizing function is underarmthrowsplineboth, called through an
% anonymous function so that the parameters are passed to it.
optfunction = @(optarg) ...
  underarmthrowsplineboth(optarg, armparms, cubicsplineparms, 0);

% Here is the optimization using fmincon (constrained optimization)
[optstar, vstar, exitflag, output] = fmincon(optfunction, optvec, [], ...
  [], [], [], lowerbound, upperbound, [], options);

figure(1); clf;
[negmaxv, bothcubic.ts, bothcubic.xs] = ...
  underarmthrowsplineboth(optstar, armparms, cubicsplineparms, 1);
% note that the results ts and xs contain time and states only for the
% knot points, so there's no need to down-sample for animation
bothcubic.maxv = -negmaxv;
% draw the animation or stop motion cartoon
figure(2); clf;
%animatemodel(xs, armparms); % uncomment if you want to see animation
cartoonmodel(bothcubic.xs, armparms);
title(sprintf('Maximal throw vopt = %4.2f', bothcubic.maxv));
bestthrows(5) = bothcubic.maxv;
timings(5) = toc;

end % double cubic splines

%% Simulation with double cubic splines and joint constraint
% Both shoulder and elbow controlled by cubic spline trajectories.
% Define a series of knot points, which are interpolated between with
% cubic spline curves. An extra constraint is added to keep angles within
% acceptable range.

if any(whichparts == 6)
tic
% Choose how many knots to include
nknots = 10; cubicsplineparms.nknotstot = 2*nknots; % two sets of knots
cubicsplineparms.splinetype = 4; % spline order is this number minus 1
cubicsplineparms.timevec = linspace(0, armparms.maxtime, nknots);

% Set up the initial guess at some knot points
shoulderknots = armparms.Tsmax * ones(1,nknots); % all maximal
elbowknots = armparms.Temax * ones(1,nknots);    % all maximal, except
% start with low force for the first part of the throw
elbowknots(cubicsplineparms.timevec < 0.8) = 0.15*armparms.Temax; 

% Perform an optimization on the spline knots to find the max velocity.
optvec = [shoulderknots elbowknots]; % this is the vector of values to be optimized

% to set up the optimization with fmincon, define lower and upper
% bounds on the optvec values:
lowerbound = optvec * 0; 
upperbound = [ones(1,nknots)*armparms.Tsmax ones(1,nknots)*armparms.Temax];

options = optimoptions(@fmincon, 'Display', 'iter', 'MaxFunEvals', 2500, ...
  'Algorithm', 'active-set');

% The optimizing function is underarmthrowsplineboth, called through an
% anonymous function so that the parameters are passed to it.
optfunction = @(optarg) ...
  underarmthrowsplineboth(optarg, armparms, cubicsplineparms, 0);

% Here is the optimization using fmincon (constrained optimization)
[optstar, vstar, exitflag, output] = fmincon(optfunction, optvec, [], ...
  [], [], [], lowerbound, upperbound, ...
  @(optarg) underarmthrowconstraint(optarg,armparms,cubicsplineparms), options);

figure(1); clf;
[negmaxv, bothconstrained.ts, bothconstrained.xs] = ...
  underarmthrowsplineboth(optstar, armparms, cubicsplineparms, 1);
% note that the results ts and xs contain time and states only for the
% knot points, so there's no need to down-sample for animation
bothconstrained.maxv = -negmaxv;
% draw the animation or stop motion cartoon
figure(2); clf;
%animatemodel(xs, armparms); % uncomment if you want to see animation
cartoonmodel(bothconstrained.xs, armparms);
title(sprintf('Maximal throw vopt = %4.2f', bothconstrained.maxv));
bestthrows(6) = bothconstrained.maxv;
timings(6) = toc;

end % double cubic splines

%% Summarize results across all methods (NaN for ones not run)
fprintf(1,'\nHere are the best throw results:\n');
fprintf(1,'%g ', bestthrows)
fprintf(1,'\nTiming of methods: ');
fprintf(1,'%g ', timings)
fprintf(1,'\n');

% Different algorithms summarized here:
% interior-point 8.66966 8.50133 8.66443 8.48988 
% sqp            8.66966 8.72297 8.73226 8.7047 8.76195 
% active-set     8.66966 8.76643 8.87026 8.7455 8.80693 7.85217

end % rununderarmthrow

%%% underarmthrowalexander (maximal joint torques)
function [negmaxv, ts, xs] = ...
  underarmthrowalexander(Tdelay, armparms, showplot)
% negmaxv = underarmthrowalexander(optvec, p, sp) returns the 
%   (negative) maximum velocity of the ball, thrown according to
%   the shoulder and elbow delay Tdelay.
%   The shoulder is maximally activated throughout, and after
%   Tdelay, the elbow is maximally activated.
% Returns negative max velocity, so that minimizing this function with
% fmincon returns a maximum throw.
%
% Input arguments are (tdelay, parms, showplot) where
%   parms contains arm model parameters,
%   and showplot = 1 (default) or 0 indicates whether
%   to plot out the simulation results automatically.
% [negmaxv, ts, xs] = ... returns the simulation time and state
%   outputs in ts, xs

if nargin < 3    % default to showing a plot unless told otherwise
  showplot = 1;
end

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; s = armparms.s;
Tsmax = armparms.Tsmax; Temax = armparms.Temax; G = armparms.G;
thsdmax = armparms.thsdmax; thedmax = armparms.thedmax;

maxtime = armparms.maxtime;

% state vector initial condition is zero
x0 = armparms.x0; % state = [qshoulder; qelbow; qdotshoulder; qdotelbow]

% simulate the whole thing
options = odeset('events', @(t,x) fevent1(t,x,armparms)); % stops when shoulder hits 180
% or ball acceleration is zero

% first simulate shoulder only, with elbow straight
[ts1,xs1] = ode45(@funderarm1, [0 Tdelay], x0,options);
%[ts1,xs1] = ode45(@(t,x) funderarm1(t,x,armparms), [0 Tdelay], x0,options);

options = odeset('events', @(t,x) fevent2(t,x,armparms)); % stops when shoulder hits 180
% or ball acceleration is zero
if max(xs1(:,1)) < pi % simulate with elbow if shoulder didn't go too far
  [ts2,xs2] = ode45(@funderarm2, [Tdelay maxtime], xs1(end,:)',options);
%  [ts2,xs2] = ode45(@(t,x) funderarm2(t,x,armparms), [Tdelay maxtime], xs1(end,:)',options);
  ts = [ts1; ts2]; xs = [xs1; xs2];
else
  ts = ts1; xs = xs1;
end

% calculate velocity from the throw
q1s = xs(:,1); q2s = xs(:,2); u1s = xs(:,3); u2s = xs(:,4);
vballx = u1s*lua.*sin(q1s) + u2s*lfa.*sin(q2s);
vbally = -u1s*lua.*cos(q1s) - u2s*lfa.*cos(q2s);
vball = sqrt(vballx.^2 + vbally.^2);
KEball = 0.5*mball*vball.*vball;
negmaxv = -max(vball); % max ball velocity

% Plot the throw
if showplot
  clf; subplot(221); plot(ts, [xs(:,1) xs(:,2)-xs(:,1)]); 
  xlabel('time'); ylabel('angles (rad)');
  subplot(222); plot(ts, [xs(:,3) xs(:,4)-xs(:,3)]); 
  xlabel('time'); ylabel('vels (rad/s)');
  subplot(223); plot(ts, vball); xlabel('time'); ylabel('ball velocity');
  hold on; plot(ts, vballx, 'g--', ts, vbally, 'r--'); 
  legend('speed','vx','vy','Location','Northwest');
  subplot(224); plot([0 ts(end)],[Tsmax Tsmax], 'b'); hold on;
  plot([0 Tdelay Tdelay ts(end)],[0 0 Temax Temax], 'g');
  xlabel('time'); ylabel('torques'); 
  legend('shoulder', 'elbow', 'Location', 'Southeast');
end

function xdot = funderarm1(t, x)
% state-derivative for under-arm throw, one-segment version
% where shoulder is maximally activated and elbow is kept at constant angle
% parameters are stored in the structure armparms

% State assignments
q1 = x(1); q2 = x(2); 
u1 = x(3); u2 = x(4); 

c1 = cos(q1); c2 = cos(q2); c1m2 = cos(q1 - q2); s1m2 = sin(q1 - q2); 

MM = zeros(2,2); rhs = zeros(2,1);

% The mass matrix is for the two-segment model, but we will lump them
% together:
MM(1,1) = mball*(lua*lua) + mfa*(lua*lua) + mua*(rua*rua); MM(1,2) = ...
c1m2*lfa*lua*mball + c1m2*lua*mfa*rfa; 
MM(2,1) = MM(1,2); MM(2,2) = mball*(lfa*lfa) + mfa*(rfa*rfa); 
k1 = rua^2*mua + 0.64*s^2*(mfa+mball);
k2 = 0.8*s*(rfa*mfa+0.8*s*mball);
k3 = rfa^2*mfa+0.64*s^2*mball;

thsd = u1;
% Calculate the torque from the torque-velocity curve for the shoulder
if thsd < thsdmax,
  Ts = Tsmax*(thsdmax-thsd)/(thsdmax+G*thsd); 
else
  Ts = 0; 
end

% since we're keeping the elbow straight, the total inertia at the shoulder
% is given here:
effM = MM(1,1) + MM(1,2) + MM(2,1) + MM(2,2); 
% This is the effective inertia for the two segmeents
% together, where u1dot = u2dot

u1dot = effM \ Ts; % simplified equation for shoulder only
u2dot = u1dot;

xdot = [x(2+1:2*2); u1dot; u2dot];

end % funderarm1

function xdot = funderarm2(t,x)
% state-derivative for under-arm throw, two-segment version
% where shoulder and elbow are both maximally activated

% State assignments
q1 = x(1); q2 = x(2); 
u1 = x(3); u2 = x(4); 

c1 = cos(q1); c2 = cos(q2); c1m2 = cos(q1 - q2); s1m2 = sin(q1 - q2); 

MM = zeros(2,2); rhs = zeros(2,1);

% Mass Matrix
MM(1,1) = mball*(lua*lua) + mfa*(lua*lua) + mua*(rua*rua); MM(1,2) = ...
c1m2*lfa*lua*mball + c1m2*lua*mfa*rfa; 
MM(2,1) = MM(1,2); MM(2,2) = mball*(lfa*lfa) + mfa*(rfa*rfa); 

thsd = u1;
% Calculate the torque from the torque-velocity curve for the shoulder
if thsd < thsdmax, 
  Ts = Tsmax*(thsdmax-thsd)/(thsdmax+G*thsd); 
else
  Ts = 0; 
end

thed = u2-u1;
% Calculate the torque from the torque-velocity curve for the elbow
if thed < thedmax, 
  Tef = Temax*(thedmax-thed)/(thedmax+G*thed); 
else
  Tef = 0; 
end

% righthand side terms
rhs(1) = c1*g*0*(lua*(mball + mfa) + mua*rua) - Tef + Ts - ...
  s1m2*lua*(lfa*mball + mfa*rfa)*(u2*u2); 
rhs(2) = c2*g*0*(lfa*mball + mfa*rfa) + Tef + s1m2*lua*(lfa*mball + ...
mfa*rfa)*(u1*u1); 
% Note that gravity is neglected, but can be turned on in these eqns
% of motion by removing the "0*" above

udot = MM\rhs;
xdot = [x(2+1:2*2); udot];

end % funderarm2

function [value, isterm, dir] = fevent1(t, x, armparms)
% stop simulation when shoulder reaches 180
% or optionally when ball acceleration goes to 0

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; s = armparms.s;

value(1) = x(1) - pi; % shoulder at 180 deg 
isterm(1) = 1;
dir(1) = 1;

value(2) = 1;

if 1 % detect ball acceleration
  xdot = funderarm1(t, x);
  u1 = xdot(1); u2 = xdot(2); u1dot = xdot(3); u2dot = xdot(4);
  c1 = cos(x(1)); s1 = sin(x(1)); c2 = cos(x(2)); s2 = sin(x(2));
  J = [lua*s1 lfa*s2; -lua*c1 -lfa*c2];
  Jdot = [lua*c1*u1 lfa*c2*u2; lua*s1*u1 lfa*s2*u2];
  aball = J*[u1dot; u2dot] + Jdot*[u1; u2];
  accball = sqrt(sum(aball.^2));
  value(2) = accball;
end

isterm(2) = 1;
dir(2) = 0;

end % fevent1

function [value, isterm, dir] = fevent2(t, x, armparms)
% stop simulation when shoulder reaches 180
% or optionally when ball acceleration goes to 0

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; s = armparms.s;

value(1) = x(1) - pi; % shoulder at 180 deg 
isterm(1) = 1;
dir(1) = 1;

value(2) = 1;
% Disabled: ball acceleration reaches 0; to enable use "if 1"
if 1 % % detect ball acceleration
  xdot = funderarm2(t, x);
  u1 = xdot(1); u2 = xdot(2); u1dot = xdot(3); u2dot = xdot(4);
  c1 = cos(x(1)); s1 = sin(x(1)); c2 = cos(x(2)); s2 = sin(x(2));
  J = [lua*s1 lfa*s2; -lua*c1 -lfa*c2];
  Jdot = [lua*c1*u1 lfa*c2*u2; lua*s1*u1 lfa*s2*u2];
  aball = J*[u1dot; u2dot] + Jdot*[u1; u2];
  accball = sqrt(sum(aball.^2));
  value(2) = accball;
end

isterm(2) = 1;
dir(2) = 0;

end % fevent2

end % underarmthrowalexander

%%% underarmthrowsplineelbow (shoulder maximal, elbow spline)
function [negmaxv, ts, xs] = ...
  underarmthrowsplineelbow(optvec, armparms, splineparms, showplot)
% negmaxv = underarmthrowsplineelbow(optvec, p, sp) returns the 
%   (negative) maximum velocity of the ball, thrown according to
%   the shoulder and elbow torque parameters specified by optvec.
% Performs a simulation based on the splines formed 
% where optvec = [elbowknots], where the shoulder is fully activated
% at all times, and the elbow torque is specified by elbowknots.
% Returns negative max velocity, so that minimizing this function with
% fmincon returns a maximum throw.
%
% Input arguments are (optvec, parms, splineparms, showplot) where
%   parms contains arm model parameters, splineparms contains info 
%   about the spline, and showplot = 1 (default) or 0 indicates whether
%   to plot out the simulation results automatically.
% [negmaxv, ts, xs] = ... returns the simulation time and state
%   outputs in ts, xs
% 

if nargin < 4    % default to showing a plot unless told otherwise
  showplot = 1;
end

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; 
Tsmax = armparms.Tsmax; Temax = armparms.Temax; G = armparms.G;
thsdmax = armparms.thsdmax; thedmax = armparms.thedmax;

maxtime = armparms.maxtime;

% state vector initial condition is zero
x0 = [0; 0; 0; 0]; % state = [qshoulder; qelbow; qdotshoulder; qdotelbow]

% Shoulder and elbow torques
% where shoulder is assumed to be maximally activated at all times
shoulderknots = Tsmax * ones(1,splineparms.nknots); 
elbowknots = optvec(1:splineparms.nknots); % and elbow defined by spline

splines.shoulder = ...
  spapi(splineparms.splinetype, splineparms.timevec, shoulderknots); % linear spline
splines.elbow = ...
  spapi(splineparms.splinetype, splineparms.timevec, elbowknots);

% simulate the whole thing
options = odeset('events', @(t,x) fevent(t,x,armparms,splines)); % stopping condition
[ts,xs] = ode45(@(t,x) funderarmspline(t,x, armparms, splines), splineparms.timevec, x0,options);

% calculate velocity from the throw
q1s = xs(:,1); q2s = xs(:,2); u1s = xs(:,3); u2s = xs(:,4);
vballx = u1s*lua.*sin(q1s) + u2s*lfa.*sin(q2s);
vbally = -u1s*lua.*cos(q1s) - u2s*lfa.*cos(q2s);
vball = sqrt(vballx.^2 + vbally.^2);
KEball = 0.5*mball*vball.*vball;
negmaxv = -max(vball); % max ball velocity

% Plot the throw
if showplot
  clf; subplot(221); plot(ts, [xs(:,1) xs(:,2)-xs(:,1)]); 
  xlabel('time'); ylabel('angles (rad)');
  subplot(222); plot(ts, [xs(:,3) xs(:,4)-xs(:,3)]); 
  xlabel('time'); ylabel('vels (rad/s)');
  subplot(223); plot(ts, vball); xlabel('time'); ylabel('ball velocity');
  hold on; plot(ts, vballx, 'g--', ts, vbally, 'r--'); 
  legend('speed','vx','vy','Location','Northwest');
  subplot(224); fnplt(splines.shoulder,'b'); hold on
  fnplt(splines.elbow,'g');
  xlabel('time'); ylabel('torques'); 
  legend('shoulder', 'elbow', 'Location', 'Southeast');
end

end %underarmthrowsplineelbow

function xdot = funderarmspline(t,x, armparms, splines)
% state-derivative for under-arm throw, two-segment version
% with spline for input.

% This is a nested function, where model parameters are defined
% in the enclosing function. These include the following:
%   mball, lua, lfa, rua, rfa, etc.

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; s = armparms.s;
Tsmax = armparms.Tsmax; Temax = armparms.Temax; G = armparms.G;
thsdmax = armparms.thsdmax; thedmax = armparms.thedmax;

maxtime = armparms.maxtime;

% Shoulder and elbow torques are the isometric torques to be applied
% and are scaled below to account for force-velocity curve.
shouldertorque = fnval(splines.shoulder,t); % splines are defined in the
elbowtorque = fnval(splines.elbow,t);       % enclosing function

% State assignments
q1 = x(1); q2 = x(2); u1 = x(3); u2 = x(4); 

c1 = cos(q1); c2 = cos(q2); c1m2 = cos(q1 - q2); s1m2 = sin(q1 - q2); 

MM = zeros(2,2); rhs = zeros(2,1);

% Mass Matrix
MM(1,1) = mball*(lua*lua) + mfa*(lua*lua) + mua*(rua*rua); MM(1,2) = ...
c1m2*lfa*lua*mball + c1m2*lua*mfa*rfa; 
MM(2,1) = MM(1,2); MM(2,2) = mball*(lfa*lfa) + mfa*(rfa*rfa); 

% Apply torque-velocity curves for shoulder and elbow
thsd = u1;
if thsd < 0            % muscle is lengthening, so set at isometric max
  Ts = shouldertorque * 1;
elseif thsd < thsdmax, % muscle is shortening at less than vmax
  Ts = shouldertorque*(thsdmax-thsd)/(thsdmax+G*thsd); 
else
  Ts = 0;              % no force beyond vmax
end

thed = u2-u1;
if thed < 0           % muscle is lengthening, so set at isometric max
  Tef = elbowtorque * 1;
elseif thed < thedmax  % muscle is shortening at less than vmax
  Tef = elbowtorque*(thedmax-thed)/(thedmax+G*thed); 
else
  Tef = 0;             % no force beyond vmax
end

% righthand side terms
rhs(1) = c1*g*0*(lua*(mball + mfa) + mua*rua) - Tef + Ts - ...
  s1m2*lua*(lfa*mball + mfa*rfa)*(u2*u2); 
rhs(2) = c2*g*0*(lfa*mball + mfa*rfa) + Tef + s1m2*lua*(lfa*mball + ...
mfa*rfa)*(u1*u1); 

udot = MM\rhs;
xdot = [x(2+1:2*2); udot];

end % funderarmlinearspline


%%% underarmthrowsplineboth (shoulder spline, elbow spline)
function [negmaxv, ts, xs] = ...
  underarmthrowsplineboth(optvec, armparms, splineparms, showplot)
% negmaxv = underarmthrowsplineboth(optvec, p, sp) returns the 
%   (negative) maximum velocity of the ball, thrown according to
%   the shoulder and elbow torque parameters specified by optvec.
% Performs a simulation based on the splines formed 
% where optvec = [elbowknots shoulderknots], where the shoulder is fully activated
% at all times, and the elbow torque is specified by elbowknots.
% Returns negative max velocity, so that minimizing this function with
% fmincon returns a maximum throw.
%
% Input arguments are (optvec, parms, splineparms, showplot) where
%   parms contains arm model parameters, splineparms contains info 
%   about the spline, and showplot = 1 (default) or 0 indicates whether
%   to plot out the simulation results automatically.
% [negmaxv, ts, xs] = ... returns the simulation time and state
%   outputs in ts, xs
% 

if nargin < 4    % default to showing a plot unless told otherwise
  showplot = 1;
end

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; 
Tsmax = armparms.Tsmax; Temax = armparms.Temax; G = armparms.G;
thsdmax = armparms.thsdmax; thedmax = armparms.thedmax;

maxtime = armparms.maxtime;

% state vector initial condition is zero
x0 = [0; 0; 0; 0]; % state = [qshoulder; qelbow; qdotshoulder; qdotelbow]

% Shoulder and elbow torques
% where shoulder is assumed to be maximally activated at all times
nknots = splineparms.nknotstot/2;
shoulderknots = optvec(1:nknots);  
elbowknots = optvec(nknots+(1:nknots)); % and elbow defined by spline

splines.shoulder = ...
  spapi(splineparms.splinetype, splineparms.timevec, shoulderknots); % linear spline
splines.elbow = ...
  spapi(splineparms.splinetype, splineparms.timevec, elbowknots);

% simulate the whole thing
options = odeset('events', @(t,x) fevent(t,x,armparms,splines)); % stopping condition
[ts,xs] = ode45(@(t,x) funderarmspline(t,x,armparms,splines),...
  splineparms.timevec, x0,options);

% calculate velocity from the throw
q1s = xs(:,1); q2s = xs(:,2); u1s = xs(:,3); u2s = xs(:,4);
vballx = u1s*lua.*sin(q1s) + u2s*lfa.*sin(q2s);
vbally = -u1s*lua.*cos(q1s) - u2s*lfa.*cos(q2s);
vball = sqrt(vballx.^2 + vbally.^2);
KEball = 0.5*mball*vball.*vball;
negmaxv = -max(vball); % max ball velocity

% Plot the throw
if showplot
  clf; subplot(221); plot(ts, [xs(:,1) xs(:,2)-xs(:,1)]); 
  xlabel('time'); ylabel('angles (rad)');
  subplot(222); plot(ts, [xs(:,3) xs(:,4)-xs(:,3)]); 
  xlabel('time'); ylabel('vels (rad/s)');
  subplot(223); plot(ts, vball); xlabel('time'); ylabel('ball velocity');
  hold on; plot(ts, vballx, 'g--', ts, vbally, 'r--'); 
  legend('speed','vx','vy','Location','Northwest');
  subplot(224); fnplt(splines.shoulder,'b'); hold on
  fnplt(splines.elbow,'g');
  xlabel('time'); ylabel('torques'); 
  legend('shoulder', 'elbow', 'Location', 'Southeast');
end

end %underarmthrowsplineboth


function [negminelbow,ceq] = ...
  underarmthrowconstraint(optvec, armparms, splineparms)
% negmaxv = underarmthrowconstraint(optvec, p, sp) returns the 
%   (negative) minimum elbow angle, thrown according to
%   the shoulder and elbow torque parameters specified by optvec.
% Performs a simulation based on the splines formed 
% where optvec = [elbowknots shoulderknots], where the shoulder is fully activated
% at all times, and the elbow torque is specified by elbowknots.
% We require that the elbow angle be non-negative and for fmincon
% the constraint is the negative of the minimum elbow angle. This is
% because fmincon maintains the constraint function as c <= 0.
%
% Input arguments are (optvec, parms, splineparms) where
%   parms contains arm model parameters, splineparms contains info 
%   about the spline.
% The output ceq is returned empty, meaning there are no equality
% constraints

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; 
Tsmax = armparms.Tsmax; Temax = armparms.Temax; G = armparms.G;
thsdmax = armparms.thsdmax; thedmax = armparms.thedmax;

maxtime = armparms.maxtime;

% state vector initial condition is zero
x0 = [0; 0; 0; 0]; % state = [qshoulder; qelbow; qdotshoulder; qdotelbow]

% Shoulder and elbow torques
% where shoulder is assumed to be maximally activated at all times
nknots = splineparms.nknotstot/2;
shoulderknots = optvec(1:nknots);  
elbowknots = optvec(nknots+(1:nknots)); % and elbow defined by spline

splines.shoulder = ...
  spapi(splineparms.splinetype, splineparms.timevec, shoulderknots); % linear spline
splines.elbow = ...
  spapi(splineparms.splinetype, splineparms.timevec, elbowknots);

% simulate the whole thing
options = odeset('events', @(t,x) fevent(t,x,armparms, splines)); % stopping condition
[ts,xs] = ode45(@(t,x) funderarmspline(t,x,armparms,splines),...
  splineparms.timevec, x0,options);

% calculate velocity from the throw
q1s = xs(:,1); q2s = xs(:,2); u1s = xs(:,3); u2s = xs(:,4);
vballx = u1s*lua.*sin(q1s) + u2s*lfa.*sin(q2s);
vbally = -u1s*lua.*cos(q1s) - u2s*lfa.*cos(q2s);
vball = sqrt(vballx.^2 + vbally.^2);

negminelbow = -min(q2s-q1s);
ceq = [];

end % constraint

function [value, isterm, dir] = fevent(t, x, armparms,splines)
% stop simulation when shoulder reaches 180
% or optionally when ball acceleration goes to 0

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; s = armparms.s;

value(1) = x(1) - pi; % shoulder at 180 deg 
isterm(1) = 1;
dir(1) = 1;

value(2) = 1;
% Disabled: ball acceleration reaches 0; to enable use "if 1"
if 1 % detect ball acceleration 0
  
  xdot = funderarmspline(t, x, armparms, splines);
  u1 = xdot(1); u2 = xdot(2); u1dot = xdot(3); u2dot = xdot(4);
  c1 = cos(x(1)); s1 = sin(x(1)); c2 = cos(x(2)); s2 = sin(x(2));
  J = [lua*s1 lfa*s2; -lua*c1 -lfa*c2];
  Jdot = [lua*c1*u1 lfa*c2*u2; lua*s1*u1 lfa*s2*u2];
  aball = J*[u1dot; u2dot] + Jdot*[u1; u2];
  accball = sqrt(sum(aball.^2));
  value(2) = accball;
end

isterm(2) = 1;
dir(2) = 0;

end % fevent


%% armenergy
function KE = armenergy(ts, xs, armparms)
% computes kinetic energy of the system, with time and states from
% a simulation arranged in arrays ts and xs, and with model parameters
% in armparms.
% Present model neglects gravity, so no potential energy is computed.

mua = armparms.mua; mball = armparms.mball; mfa = armparms.mfa;
lua = armparms.lua; lfa = armparms.lfa;
rua = armparms.rua; rfa = armparms.rfa;
g = armparms.g; s = armparms.s;
Tsmax = armparms.Tsmax; Temax = armparms.Temax; G = armparms.G;
thsdmax = armparms.thsdmax; thedmax = armparms.thedmax;

KE = ts * 0; % allocate space for a vector of kinetic energies
for i = 1:length(ts)

  % State assignments
  q1 = xs(i,1); q2 = xs(i,2);
  u1 = xs(i,3); u2 = xs(i,4);
  
  c1 = cos(q1); c2 = cos(q2); c1m2 = cos(q1 - q2); s1m2 = sin(q1 - q2);
  
  MM = zeros(2,2); rhs = zeros(2,1);
  
  % The mass matrix is for the two-segment model, but we will lump them
  % together:
  MM(1,1) = mball*(lua*lua) + mfa*(lua*lua) + mua*(rua*rua); MM(1,2) = ...
    c1m2*lfa*lua*mball + c1m2*lua*mfa*rfa;
  MM(2,1) = MM(1,2); MM(2,2) = mball*(lfa*lfa) + mfa*(rfa*rfa);
  
  KE(i) = 0.5*[u1 u2]*MM*[u1; u2];

end

end % armenergy

%% Helper functions for drawing and animation
function hfigure = drawmodel(x, hfigure, posn, armparms)
% draws the arm with state x, using figure handle hfigure,
% at position posn, with parameters armparms

lua = armparms.lua; lfa = armparms.lfa;

ths = x(1); thf = x(2); % shoulder and forearm angles

armx = [0 -lua*cos(x(1)) -lua*cos(x(1))-lfa*cos(x(2))] + posn(1);
army = [0 -lua*sin(x(1)) -lua*sin(x(1))-lfa*sin(x(2))] + posn(2);

if isempty(hfigure) % no figure handle given
  hfigure = plot(armx, army, 'b-o'); 
else
  set(hfigure, 'XData', armx, 'YData', army);
end


end % drawmodel

function cartoonmodel(xs, armparms)
% draw the model across the page to show a cartoon of its motion

dx = 1.1; % advance model horizontally by this much each frame

clf; hold on; axis equal

for i = 1:size(xs,1)
  xposn = (i-1)*dx;
  drawmodel(xs(i,:), [], [xposn 0], armparms);
end

end % cartoonmodel

function animatemodel(xs, armparms)
% animate the model 

clf; hold on; axis equal; set(gca,'xlim',[-2.1 2.1],'ylim',[-2.1 2.1]);   
h = [];

for i = 1:size(xs,1)
  h = drawmodel(xs(i,:), h, [0 0], armparms);
  pause(0.1)
end

end % animatemodel