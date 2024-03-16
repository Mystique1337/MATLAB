clc;
clear all
close all

tic
%% Code to predict COVID-19 transmission in a classroom

% For additional information on some of these estimates see https://tinyurl.com/covid-estimator

%% constants
rho = 1000; % density of water, kg/m3 (assumed density of exhaled breath droplets)
g = 9.81; % gravitational acceleration, m/s2
mu = 1.8e-5; % viscosity of air, kg/m-s
lambda = 65e-9; % mean free path of air, m
kB = 1.38e-23; % Boltzmann constant, J/K
R = 8.314; % ideal gas constant, J/mol-K

%% environmental parameters
T = 72; % temperature of room, F
T = (5/9)*(T-32); % temperature of room, C
RH = 50; % indoor RH, assuming the upper bound of indoor recommendations
S = 0.05; % solar irradiation, W/m2, assuming a small value (0.05?) indoors
decay = k_inf(T,RH,S)/60; % 1st-order decay constant in air, 1/s 
% See https://www.tandfonline.com/doi/full/10.1080/02786826.2020.1829536


%% facility parameters
A = 100; % assumed floor area, m2 
h = 3; % assumed height, m
SA = 2*A + 4*sqrt(A)*h; % calculated surface area (walls, ceiling, and floor), m2
V = A*h; % volume of room, m3
ACH = 5/3600; % air changes per hour, converted to per second (3/hr is "typical" for a university classroom)

%% human behavior
Q = 16/24/3600; % breathing rate, m3/s (assuming 16 m3/day)
people = 46; % number of people in the room (default is 50, the maximum occupancy)
mask_out_eff = 0.50; % everyone is wearing a "standard" mask that removes 50% of the exhaled droplets
mask_in_eff = 0.30; % everyone is wearing a "standard" mask that removes 30% of the inhaled droplets
duration = 55; % duration of one class meeting
classes = 42; % number of in-person classes during AU21

% The term "quanta" is used when dose-response is unknown. This is what
% people have been using for COVID-19. One recent estimate is that 18.6
% quanta are emitted per hour for the OG variant, but this is doubled for
% the Delta variant.
q = 2*18.6;

%% treatment parameters
% some will be removed by central HVAC
filter_eff = 0.90; % assume recirculated air filtered at 90% efficiency with 
% "MERV-13" filter

% we could remove more if we used a portable air cleaner within the room
nC = 0; % the number of portable air cleaners
CADR = 800*0.000472; % clean air delivery rate (CADR), m3/s (0.000472 converts cfm to m3/s)
% See https://www.epa.gov/indoor-air-quality-iaq/air-cleaners-and-air-filters-home
% for guidance on selection -- dependent on floor area in room
% We have two of these in our classroom, and I've conservatively used CADR
% = 800
CADR_eff = 0.99; % assume 99% efficiency because these have HEPA filters

% Some have proposed UV disinfection of air. how this works could be
% uncertain. We can ignore (?) since this is built into the new decay
% calculator.
UV = 0/60; % UV disinfection rate, per second (0.2/min default?) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Everything below here is the "guts" of the model. %%%%%%%%%%%%%%%
%%%%%%%%% If you change something there, you could "break" it.%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% aerosol parameters
dp_min = 10e-9; % minimum particle diameter, m
mp_min = rho*pi/6*(dp_min)^3; % minimum particle mass, kg

bins = 50; % this defines 50 particle size bins
scale = 2; % this is the factor by which mass increases as we go from bin n to n+1
mp_max = mp_min*scale^bins; % maximum particle mass, kg
dp_max = (6*mp_max/rho/pi)^(1/3); % maximum particle diameter, m

% The following are size distribution parameters, to match with Johnson et al.  
% for breathing (see https://doi.org/10.1016/j.jaerosci.2011.07.009)
% this is a tri-modal distribution with dpgi as the geometric mean of mode
% i, sgi is the geometric standard deviation of mode i, and Ni is the
% particle number concentration in mode i

dpg1 = 1.6e-6; % m
sg1 = 1.30;
N1 = 0.069; % #/cm3
dpg2 = 2.6e-6;
sg2 = 1.66;
N2 = 0.085; % #/cm3
dpg3 = 145e-6;
sg3 = 1.795;
N3 = 0.001;

% this loop takes everyhing and creates the particle size (and mass) distributions
for i = 1:bins+1
    mp(1) = mp_min;
    mp(i) = mp_min*scale^(i-1);
    
    dp0(1) = dp_min;
    dp0(i) = (6*mp(i)/rho/pi)^(1/3);
    
    N0(i) = (N1/(sqrt(2*pi)*log10(sg1)))*exp(-(log10(dp0(i)) - log10(dpg1))^2/(2*log10(sg1)^2))...
    + (N2/(sqrt(2*pi)*log10(sg2)))*exp(-(log10(dp0(i)) - log10(dpg2))^2/(2*log10(sg2)^2))...
    + (N3/(sqrt(2*pi)*log10(sg3)))*exp(-(log10(dp0(i)) - log10(dpg3))^2/(2*log10(sg3)^2));

    N01(i) = (N1/(sqrt(2*pi)*log10(sg1)))*exp(-(log10(dp0(i)) - log10(dpg1))^2/(2*log10(sg1)^2));
    N02(i) = (N2/(sqrt(2*pi)*log10(sg2)))*exp(-(log10(dp0(i)) - log10(dpg2))^2/(2*log10(sg2)^2));
    N03(i) = (N3/(sqrt(2*pi)*log10(sg3)))*exp(-(log10(dp0(i)) - log10(dpg3))^2/(2*log10(sg3)^2));

    M0(i) = N0(i)*1e6*rho*1e9*pi/6*dp0(i)^3;
end

N0_int = trapz(dp0,N0*1e6); % number concentration of exhaled droplets (#/m3)

% figure
% loglog(dp0*1e6,N01,'r','LineWidth',2)
% hold on
% loglog(dp0*1e6,N02,'b','LineWidth',2)
% loglog(dp0*1e6,N03,'g','LineWidth',2)
% loglog(dp0*1e6,N0,'k--','LineWidth',3)
% set(gca,'FontSize',20)
% xlabel('Diameter (\mum)','FontSize',30)
% ylabel('dN/dlog_1_0(d_p) (# cm^-^3 \mum^-^1)','FontSize',30)
% axis([1e-2 1e3 1e-5 1e0])
% legend('Mode 1','Mode 2','Mode 3','Total','Location','NorthWest')

%% Implementing the box model

% we first need to define some emission rates
E0 = N0*1e6*Q; % droplet emission rate for breathing, #/s
scale = 0.5e4; % scaling up breathing emissions to speaking, based on 
    % https://doi.org/10.1371/journal.pone.0227699

% total emission rate, assuming that one person is talking at a time, #/s
E = E0*(people-1)*(1-mask_out_eff) + N0*1e6*Q*scale*(1-mask_out_eff); 



dt = 60; % time step for model output, seconds
time = 0:dt:60*180; % time array, from 0 to duration of class, seconds

for j = 2:length(time) % do calculations spanning the 2nd through final entry
    % in the time array
    for i = 1:bins+1 % do calculations across all particle sizes
        N(i,1) = 0*N0(i); % initial number of droplets associated with each 
        % size, assumed to be 0. We can think of this as a "number
        % distribution" but it's units are #, NOT #/cm3
        
        kw(i) = particle_loss(dp0(i),ACH,V,T); % particle wall loss rate, 1/s
        % This is calculated using a function at the bottom of this code
        % and is, in part, a function of settling velocities and
        % diffusivities. Hence, it is a function of particle size, so we do
        % this inside the loop.
        
        % This is a "stiff" system, and I wrote this code before I
        % re-taught myself a Matlab function that can handle this type of
        % ODE. Instead, this uses "sub-steps" to avoid the need for using 1
        % ms time steps above (e.g., there would be 4.8 million outputs
        % per size bin because there are 4.8 million milliseconds in 80 
        % minutes.
        substep = 1e-3; % seconds
        time_step = time(j-1); % this is a time step that gets over-written 
        % during each sub-step
        
        N_old(i) = N(i,j-1); % this is a "dummy" variable for the number 
        % distribution that gets over-written during each sub-step

        
        % This while loops basically says "Stay in this loop until you have
        % done dt/substep calculations."
        while time(j) > time_step
            % Here is where we actually do the box model calculations. This
            % uses a backwards difference scheme (see PowerPoint file).
            % Like N_old, N_new is a number distribution (#)
            N_new(i) = N_old(i) + substep*(E(i) - N_old(i)*kw(i) - N_old(i)*ACH*filter_eff - N_old(i)*decay - N_old(i)*nC*CADR/V*CADR_eff - N_old(i)*UV);
            % Go forward one sub-step
            time_step = time_step + substep;
            % Update the number distribution
            N_old(i) = N_new(i);
        end
        
        N(i,j) = N_old(i); % stored number distribution for time step j (N)
        
    end
end


%% Making figures

% Note: see Matlab help for explanation on plotting functions. You can
% type, e.g. "doc loglog", for the help file on the loglog plot function.

% This figure plots size-resolved sources and sinks
figure
loglog(dp0*1e6,E*3600*.1,'c','LineWidth',3) % emissions
hold on
loglog(dp0*1e6,(kw+ACH+decay+CADR/V*CADR_eff+UV)*3600,'r','LineWidth',3) % total removal
loglog(dp0*1e6,(kw)*3600,'k--','LineWidth',3) % wall loss only
loglog(dp0*1e6,ACH*ones(size(dp0))*3600,'b--','LineWidth',3) % air exchange only
loglog(dp0*1e6,decay*ones(size(dp0))*3600,'g--','LineWidth',3) % airborne decay
set(gca,'FontSize',20)
xlabel('Diameter (\mum)','FontSize',30)
ylabel('Rate (# hr^-^1)','FontSize',30)
axis([1e-2 1e3 1e-3 1e8])
legend('Emission (+)','Total loss (-)','Particle motion (-)','Air exchange (-)','Decay (-)','Location','NorthWest')


% This figure is a contour plot that provides a time series of the particle
% size distribution within the room
N(:,1) = 0*N0'; % We don't actually have an entry of N for the first time
% step. This takes care of that.
figure
contourf(time/60,dp0*1e6,(N/V),20,'LineStyle','none')
set(gca,'FontSize',20)
set(gca,'YScale','log')
colorbar EastOutside
axis([0 time(end)/60 dp_min*1e6 dp_max*1e6])
xlabel('Elapsed room (min)','FontSize',30)
ylabel('Particle diameter (\mum)','FontSize',30)
cbar = colorbar;
set(get(cbar,'label'),'String','dN/dlog(d_p) (# m^-^3 \mum)','FontSize',24)

% We need a few intermediate calculations for the next two figures
q = q/3600; % "quanta" emitted per second in room
qN = q/trapz(dp0,E); % "quanta" per particle
number = trapz(dp0,N)/V; % total number concentration of exhaled droplets (#/m3)
quanta = qN*number; % total quanta concentration (#/m3)

% Plotting the number concentration of exhaled droplets
figure
plot(time/60,number,'LineWidth',3)
set(gca,'FontSize',20)
xlabel('Time (min)','FontSize',30)
ylabel({'Droplet concentration'; 'in air (# m^-^3)'},'FontSize',30)

% Plotting the quanta concentration
figure
plot(time/60,quanta,'LineWidth',3)
set(gca,'FontSize',20)
xlabel('Time (min)','FontSize',30)
ylabel({'Quanta concentration'; 'in air (# m^-^3)'},'FontSize',30)

% We need some more intermediate calculations to calculate risk
V_inh = Q*time; % total volume inhaled at a given time step
qinh = V_inh.*quanta*(1-mask_in_eff); % total quanta inhaled at a given time step
risk = 1 - exp(-qinh); % risk of infection estimate based on Wells-Riley model

% Plotting the probability of a new case
figure
plot(time/60,risk,'LineWidth',3)
set(gca,'FontSize',20)
xlabel('Time (min)','FontSize',30)
ylabel({'Probability of new case'},'FontSize',30)

% Scaling the probability by number of people to estimate the number of new
% COVID-19 cases per class session
figure
plot(time/60,risk*people,'LineWidth',3)
set(gca,'FontSize',20)
xlabel('Time (min)','FontSize',30)
ylabel({'Number of new'; 'cases per class'},'FontSize',30)

new_cases = risk*(people-1);
new_cases_per_class = new_cases(duration+1) % duration + 1 pulls the value associated with the class duration

% semester = 1 - (1-risk(56))^classes % probability of infection at some point during the semester

toc
%% Viral decay rate function
function decay = k_inf(T,RH,S)
decay = 0.1603 + 0.04018*((T-20.615)/10.585) + 0.02176*((RH-45.235)/28.665) ...
    + 0.14369*((S-0.95)/0.95) + 0.02636*((T-20.615)/10.585)*((S-0.95)/0.95);
end

%% Particle loss coefficient function
function kw = particle_loss(dp,ACH,V,T)

rho = 1000; % kg/m3
g = 9.81;
mu = 1.8e-5;
lambda = 65e-9;
kB = 1.38e-23;
T = T + 273;

ke = 0.6*ACH + .25/3600; % turbulent kinetic energy in room, based on https://dx.doi.org/10.1021/es103080p


L = V^(1/3); % length scale for loss

Cc = 1 + 2*lambda/dp*(1.257 + 0.4*exp(-1.1*dp/(2*lambda)));
vts = rho*dp^2*g*Cc/(18*mu);
D = kB*T*Cc/(3*pi*mu*dp);
x = pi*vts/(2*sqrt(ke*D));

kw = 1/(L)*(8*sqrt(ke*D)/pi + vts*coth(x/2)); % following https://doi.org/10.1080/02786828308958636

end

