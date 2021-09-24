%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of ODEs describing enzymatic reaction kinetics in the presence of a competitive inhibitor:
% E + S <-> ES      (k1f, k1r)
% ES -> E + P       (k2f)
% E + I <-> EI      (kif, kir) 
% EI + S <-> ESI    (k1f, k1r)
% ES + I <-> ESI    (kif, kir)

clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define parameters

k1f = 0.01;         % [uM^-1 s^-1]
k1r = 1;            % [s^-1]
k2f = 0.5;          % [uM^-1 s^-1]
kif = 0.01;         % [uM^-1 s^-1]
kir = 1;            % [s^-1]

P0 = 0;     % [uM] initial concentration of product
ES0 = 0;
Etot = 10;   % [uM] total (initial) concentration of enzyme
I0 = 100;     % [uM] initial concentration of inhibitor
EI0 = 0;
S0 = 100;     % [uM] initial concentration of substrate
EIS0 = 0; 

parameters = {k1f, k1r, k2f, kif, kir};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run single simulation
y0 = [  % initial conditions of the experiment/simulation
    P0
    ES0
    Etot
    I0
    EI0
    S0
    EIS0
    ];
tspan = [0 5];     % [s] timespan of experiment/simulation
options = [];
[t,y] = ode15s(@noncompetitive_inh_enzyme_react_ODEfun,tspan,y0,options,parameters);

clear reaction_rate;

reaction_rate = k2f*(y(:,2)); % reaction rate = v = d[P]/dt = k2f*(ES)

% plot time-course concentrations and reaction rate
figure(1); hold on;

% concentrations
subplot(2,1,1);
plot(t,y(:,1:7),'linewidth',2);
hold on;
% set(gca,'Fontsize',15);
xlabel('Time (s)');
ylabel('Concentration (t) [\muM]');
legend('[P]','[ES]','[E]','[I]','[EI]','[S]','[EIS]');
title(['S0 = ' char(num2str(S0)) ' \muM, Etot = ' char(num2str(Etot)) ' \muM, I0 = ' char(num2str(I0)) ' \muM']);

% reaction rate
subplot(2,1,2);
plot(t,reaction_rate,'linewidth',2);
hold on;
% set(gca,'Fontsize',15);
xlabel('Time (s)');
ylabel('reaction rate (t) (\muM/s)');
legend('\nu(t) = d[P]/dt');
title(['S0 = ' char(num2str(S0)) ' \muM, Etot = ' char(num2str(Etot)) ' \muM, I0 = ' char(num2str(I0)) ' \muM']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulations for initial rate experiments --> reaction rate (v) versus initial substrate concentration
% assume (t = 1 s) is an appropriate time for the initial rate experiment, where quasi-steady-state assumptions holds

clear h1 h2 inhibitor_conditions;

I0_range = [0 30 120 480 1200];
S0_range = [120 240 480 1200 3000 10000];

for j = 1:length(I0_range)
    
    clear initial_rate;
    for i = 1:length(S0_range)
        S0 = S0_range(i);
        I0 = I0_range(j);
        parameters = {k1f, k1r, k2f, kif, kir};
        y0 = [  % initial conditions of the experiment/simulation
            P0
            ES0
            Etot
            I0
            EI0
            S0
            EIS0
            ];
        tspan = [0 5];
        options = [];
        [t,y] = ode15s(@noncompetitive_inh_enzyme_react_ODEfun,tspan,y0,options,parameters);
        
        clear reaction_rate time_steps temp;
        reaction_rate = k2f*(y(:,2)); % reaction rate = v = d[P]/dt = k2f*(ES)
        %     time_steps = 1:length(t);
        %     temp = time_steps(t > 1);
        %     initial_rate(i) = reaction_rate(temp(1)); % reaction rate at ~1 s for each substrate concentration
        initial_rate(i) = reaction_rate(end);
    end
    
    figure(2);
    hold on;
    h1(j) = plot(S0_range,initial_rate,'linewidth',2,'marker','o');
%     set(gca,'Fontsize',15);
    xlabel('[S] (\muM)');
    ylabel('\nu (\muM/s)');
    
    
    figure(3);
    hold on;
    h2(j) = plot(1./S0_range,1./initial_rate,'linewidth',2,'marker','o');
%     set(gca,'Fontsize',15);
    xlabel('1/[S] (1/\muM)');
    ylabel('1/\nu (s/\muM)');
    title('Lineweaver-Burk plot');
    
    inhibitor_conditions{j} = ['noncompetetive inh, ' char(num2str(I0_range(j))) ' \muM'];
    
end
legend(h1(:),inhibitor_conditions);
legend(h2(:),inhibitor_conditions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inhibitor dose-response relationship and potency (in the presence of a fixed concentration of substarte)

S0_range = [150 1500 15000];    % [uM] substrate concentration
I0_range = [0 10.^(-2:1:6)];    % [uM] inhibitor concentration

clear h3 h4 substrate_conditions;

for i = 1:length(S0_range)
    clear initial_rate product_conc;
    for j = 1:length(I0_range)
        I0 = I0_range(j);
        S0 = S0_range(i);
        parameters = {k1f, k1r, k2f, kif, kir};
        y0 = [  % initial conditions of the experiment/simulation
            P0
            ES0
            Etot
            I0
            EI0
            S0
            EIS0
            ];
        tspan = [0 5];
        options = [];
        [t,y] = ode15s(@noncompetitive_inh_enzyme_react_ODEfun,tspan,y0,options,parameters);
        
        clear reaction_rate time_steps temp;
        reaction_rate = k2f*(y(:,2)); % reaction rate = v = d[P]/dt = k2f*(ES)
        %     time_steps = 1:length(t);
        %     temp = time_steps(t > 1);
        %     initial_rate(i) = reaction_rate(temp(1)); % reaction rate at ~1 s for each substrate concentration
        initial_rate(j) = reaction_rate(end);
        product_conc(j) = y(end,1);
    end
    
    figure(4);
    subplot(1,3,1);
    hold on;
    % normalized reaction rate as a function of inhibitor concentration
    h3(i) = plot(log10(I0_range),initial_rate/initial_rate(1),'linewidth',2,'marker','o');
%     set(gca,'Fontsize',15);
    xlabel('Log_1_0([I])');
    ylabel('Normalized rate: \nu(I) / \nu(I=0)');
    h6 = plot([-2 6],[0.5 0.5],'--k','linewidth',2);
    title('dose-response relationship');
    
    
    
    subplot(1,3,2);
    hold on;
    % normalized product concentration (at t = tspan) as a function of inhibitor concentration
    h4(i) = plot(log10(I0_range),product_conc,'linewidth',2,'marker','o');
%     set(gca,'Fontsize',15);
    xlabel('Log_1_0([I])');
    ylabel('Product Concentration [\muM]');
    title('dose-response relationship');
    
    subplot(1,3,3);
    hold on;
    % normalized product concentration (at t = tspan) as a function of inhibitor concentration
    h5(i) = plot(log10(I0_range),product_conc/product_conc(1),'linewidth',2,'marker','o');
%     set(gca,'Fontsize',15);
    xlabel('Log_1_0([I])');
    ylabel('Normalized Product Conc.');
    h6 = plot([-2 6],[0.5 0.5],'--k','linewidth',2);
    title('dose-response relationship');
    
    substrate_conditions{i} = ['Substrate, ' char(num2str(S0_range(i))) ' \muM'];
    
end
legend(h3(:),substrate_conditions);
legend(h4(:),substrate_conditions);
legend(h5(:),substrate_conditions);






