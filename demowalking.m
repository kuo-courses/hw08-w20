%% Demonstration of class library for walking simulations
% Matlab's object oriented programming system is used to organize
% a consistent data structure for simulations, operated on by
% common functions (called methods). This file demonstrates use of
% the class library as well as the publish command.
%
% Make sure the walkclass directory is in your path.

%% Simplest walking model 2-D: walksw2
% First construct a gait with some default parameters.
w = walksw2; 
x0 = [0.0307 2*0.0307 -0.03256 -6.13654e-005];
% (the default P is very small, and you need a very good initial guess)
[w,cnvrg] = gradsearch(w, x0); % gradient method search for fixed point
onestep(w); % this is how to view the angles in a figure

%%% This is a normal gait
% Stored gait for easy retrieval, speed 1.25 m/s, 1.8 Hz
w = walksw2('normal'); % speed 1.25 m/s, 1.8 Hz
% here's how to display info about the gait:
w
% and here's how to retrieve specific info:
x0 = get(w, 'xstar');
% and here's how to do a step using a different initial condition
% (if you leave off the second argument, it defaults to the 
% fixed point xstar)
x0(1) = x0(1)*1.1; % alter one state to be 10% bigger
[xe,te,x,t] = onestep(w, x0); % xe, te are end state and time
% and x, t contain the entire trajectories.  
% since we are not starting at the fixed point, the step
% does not reproduce the initial conditions:
onestep(w, x0) - x0

%%% A slow gait
% Another gait that's already stored is a slow one, with no hip spring:
w = walksw2('slow'); % dimensionless speed 0.1, Kp = 0
[v, sl, sf] = gaitspeed(w); % this retrieves the speed, step length, step frequency
% suppose we wish to find the normal gait above.  that is done by setting
% target speed and step length, and using findgaitspeed:
v = 0.4; sf = 1.8/sqrt(9.81); sl = v/sf;
% start with gait w, and find a gait wnew that satisfies desired speed,
% step length:
[wnew, cnvrg] = findgaitspeed(w, v, sl); 
wnew

%% Anthropomorphic walking model 2D: walk2
% This is a model with human-like mass properties for the legs.
% It will be demonstrated to perform parameter studies.

w = walk2; w = gradsearch(w);
v = 0.4; sf = 1.8/sqrt(9.81); sl = v/sf;
[wnew,cnvrg] = findgaitspeed(w, v, sl);

%% Parameter study 1: Changing Kp
% What happens if the hip spring changes?
w = walk2('slow');
studyKp.Kps = [0 logspace(-3,0,10)]; 
% Perform a parameter study using specified values of parameter
studyKp.ws = parmstudy1d(w, studyKp.Kps, 'Kp');

%%% Examine parameter study results for Kp
% Post-process to get speeds and step lengths
clf; hold on;
for i = 1:length(studyKp.Kps)
  [studyKp.speeds(i), studyKp.sls(i), studyKp.sfs(i)] = ...
    walkspeed(studyKp.ws(i));
  onestep(studyKp.ws(i)); % with no output arguments, plot the angles
  [xc,tc,xs,ts,energies] = onestep(studyKp.ws(i)); % get energy info
  studyKp.energies(i) = energies;
end
snapnow; % include the previous plot in the published file

%% Parameter study 2: Changing push-off P
w = walk2('slow');
studyP.Ps = [get(w,'P') logspace(log10(0.02), log10(0.4), 15)];
studyP.ws = parmstudy1d(w, studyP.Ps, 'P');

%%% Examine parameter study rsults for P
clf; hold on;
for i = 1:length(studyP.Ps)
  [studyP.speeds(i), studyP.sls(i), studyP.sfs(i)] = ...
    walkspeed(studyP.ws(i));
  onestep(studyP.ws(i)); % with no output arguments, plot the angles
  [xc,tc,xs,ts,energies]=onestep(studyP.ws(i)); % get energy info
  studyP.energies(i) = energies;
end
snapnow; % include the previous plot in the published file

%% Summarize parameter study effects on speed, step length & frequency
clf;
subplot(231)
plot(studyKp.Kps, studyKp.speeds)
xlabel('Kp'); ylabel('Speed'); title('Varying Kp');
subplot(232)
plot(studyKp.Kps, studyKp.sls)
xlabel('Kp'); ylabel('Step length'); title('Varying Kp');
subplot(233)
plot(studyKp.Kps, studyKp.sfs)
xlabel('Kp'); ylabel('Step frequency'); title('Varying Kp');
subplot(234)
plot(studyP.Ps, studyP.speeds)
xlabel('P'); ylabel('Speed'); title('Varying P');
subplot(235)
plot(studyP.Ps, studyP.sls)
xlabel('P'); ylabel('Step length'); title('Varying P');
subplot(236)
plot(studyP.Ps, studyP.sfs)
xlabel('P'); ylabel('Step frequency'); title('Varying P');

%% Summarize parameter study effects on collision losses
% This shows how push-off work affects speed when Kp or P is
% being varied.
clf;
subplot(121)
plot([studyKp.energies.pushoffwork], studyKp.speeds); 
xlabel('Push-off work'); ylabel('Speed'); title('Varying Kp');
subplot(122)
plot([studyP.energies.pushoffwork], studyP.speeds); 
xlabel('Push-off work'); ylabel('Speed'); title('Varying P');

