clear all;
clc;

%% Problem 1 

% Part a

CPA_conc_init = [0.3 0.5 0.8 1.5 3 5 7]; % substrate 
P450_2B1_rat = [5.82 9.03 12.7 17.1 20.2 27.8 31.5]; 
P450_2B1_variant = [17.5 24.5 24.0 23.9 27.3 33.1 27.7];

% Best fit line
data_rat = polyfit(CPA_conc_init, CPA_conc_init./P450_2B1_rat,1);
Yintercept_rat = data_rat(2); 
Xintercept_rat = -data_rat(2) / data_rat(1);

data_variant = polyfit(CPA_conc_init, CPA_conc_init./P450_2B1_variant,1);
Yintercept_variant = data_variant(2); 
Xintercept_variant = -data_variant(2) / data_variant(1);

substrate = [-3:1:8]; 
fit_rat = data_rat(1)*substrate + data_rat(2);
fit_variant = data_variant(1)*substrate + data_variant(2);

Vmax_rat = 1/data_rat(1);
Vmax_variant = 1/data_variant(1);


figure(1) 
plot(CPA_conc_init, CPA_conc_init./P450_2B1_rat,'o','LineWidth',2,'Color',[0.0 0.5 0.5])
hold on
plot(CPA_conc_init,CPA_conc_init./P450_2B1_variant,'o','LineWidth',2,'Color',[0.5 0.0 0.5])
plot(substrate, fit_rat,'--','LineWidth',2,'Color',[0.0 0.8 0.8]) 
plot(substrate,fit_variant,'--','LineWidth',2,'Color',[0.8 0.0 0.8])
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
xlabel('CPA initial concentration [S] (mM)')
ylabel('Velocity of reaction / CPA initial concentration [S]/[v] (min*lit/1000)')
title('Hanes-Woolf Plot')
legend('P450 2B1 Rat Experimental Data','P450 2B1 Variant Experimental Data',...
    'P450 2B1 Rat Fitted Line','P450 2B1 Variant Fitted Line')
hold off 

%% Part b 
clear substrate fit_rat fit_variant
substrate = [0:1:7]; 
fit_rat = data_rat(1)*substrate + data_rat(2);
fit_variant = data_variant(1)*substrate + data_variant(2);

figure(2) 
hold on
plot(substrate, fit_rat,'-','LineWidth',2,'Color',[0.0 0.8 0.8]) 
plot(substrate,fit_variant,'-','LineWidth',2,'Color',[0.8 0.0 0.8])
xline(0.1,'--r');  %limit 1
xline(0.2,'--r');  %limit 2
xlabel('CPA initial concentration [S] (mM)')
ylabel('Velocity of reaction / CPA initial concentration [S]/[v] (min*lit/1000)')
title('Hanes-Woolf Plot')
legend('P450 2B1 Rat Experimental Data','P450 2B1 Variant Experimental Data',...
    '[S]=0.1 mM','[S]=0.2 mM')
hold off 
set(gca,'FontSize',14)

%% Problem 2

%% Part a 

% Lineweaver Burk Plot
substrate = [0.25 0.40 0.50 0.60 0.75 1.00];

% Row 1: I=0mM ; Row 2: I=1mM ; Row 3: I=2mM
rate_reaction = [1.00 1.70 1.90 2.10 2.40 2.50
    0.65 1.10 1.30 1.40 1.80 2.20 
    0.55 0.90 1.00 1.30 1.40 1.80];

ones(6);
figure(3)
hold on;
for i = 1:1:3
    plot(ones./substrate, ones./rate_reaction(i,:),'Linewidth',2,'Color', i*[0.0 0.3 0.3])
end
hold off;

xlabel('1/[S]')
ylabel('1/v')
legend('I = 0mM','I = 1 mM', 'I = 2 mM','Location','northwest')
title('Lineweaver-Burk Plot')
set(gca,'FontSize',14)


% Best fit line
ones(6);
data_I_0 = polyfit(ones./substrate, ones./rate_reaction(1,:),1);
data_I_1 = polyfit(ones./substrate, ones./rate_reaction(2,:),1);
data_I_2 = polyfit(ones./substrate, ones./rate_reaction(3,:),1);

Yintercept_I_0 = data_I_0(2); 
Xintercept_I_0 = -data_I_0(2) / data_I_0(1);


Yintercept_I_1 = data_I_1(2); 
Xintercept_I_1 = -data_I_1(2) / data_I_1(1);

Yintercept_I_2 = data_I_2(2); 
Xintercept_I_2 = -data_I_2(2) / data_I_2(1);

clear substrate
substrate = [-0.5 0.2]; 
one_substrate = ones(1,size(substrate,2))./substrate;
fit_I_0 = data_I_0(1)*one_substrate + data_I_0(2);
fit_I_1 = data_I_1(1)*one_substrate + data_I_1(2);
fit_I_2 = data_I_2(1)*one_substrate + data_I_2(2);

ones(size(substrate,2)); 
substrate_1 = [0.25 0.40 0.50 0.60 0.75 1.00];

figure(2) 
hold on;
plot(ones./substrate_1, ones./rate_reaction(1,:),'o','Linewidth',2,'Color', 1*[0.0 0.3 0.3])
plot(ones./substrate_1, ones./rate_reaction(2,:),'o','Linewidth',2,'Color', 2*[0.0 0.3 0.3])
plot(ones./substrate_1, ones./rate_reaction(3,:),'o','Linewidth',2,'Color', 3*[0.0 0.3 0.3])
plot(one_substrate, fit_I_0,'--','LineWidth',2,'Color',[0.0 0.3 0.3]) 
plot(one_substrate, fit_I_1,'--','LineWidth',2,'Color',2*[0.0 0.3 0.3]) 
plot(one_substrate, fit_I_2,'--','LineWidth',2,'Color',3*[0.0 0.3 0.3]) 
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
hold off;
xlabel('1/[S] (1/mM)')
ylabel('1/v (l.h/mM)')
xlim([-2,5])
legend('I = 0mM','I = 1 mM', 'I = 2 mM','I = 0mM Fitted','I = 1 mM Fitted'...
    ,'I = 2 mM Fitted','Location','northwest')
title('Lineweaver-Burk Plot')
set(gca,'FontSize',14)

%% Problem 3 

clear substrate rate_reaction ones_substrate ones_reaction
substrate = [10 20 30 50 60 80 90 110 130 140 150];
rate_reaction = [5 7.5 10 12.5 13.7 15 15 12.5 9.5 7.5 5.7];

ones_substrate = ones(1,size(substrate,2))./substrate;
ones_rate_reaction = ones(1,size(substrate,2))./rate_reaction;
rate_substrate = rate_reaction./substrate;
substrate_rate = substrate./rate_reaction;

% Best fit plot 

substrate_fit = [10 20 30 50 60 80];
rate_reaction_fit = [5 7.5 10 12.5 13.7 15];

ones_substrate_fit = ones(1,size(substrate_fit,2))./substrate_fit;
ones_rate_reaction_fit = ones(1,size(substrate_fit,2))./rate_reaction_fit;

fit = polyfit(ones_substrate_fit, ones_rate_reaction_fit(1,:),1);
data_fit = fit(1)*ones_substrate + fit(2);

Yintercept_fit = fit(2); 
Xintercept_fit = -fit(2) / fit(1);

Vmax_fit = 1/Yintercept_fit;
Km_fit = -1/Xintercept_fit;

% Lineweaver-Burk Plot
figure(4)
plot(ones_substrate,ones_rate_reaction,'-o','LineWidth',2)
hold on 
plot(ones_substrate, data_fit,'LineWidth',2)
xlabel('1/[S] (l/mg)','FontSize',14)
ylabel('1/v (l.h/mg)','FontSize',14)
title('Lineweaver-Burk','LineWidth',2,'Color','k','FontSize',14)

% Calculating KSI
Smax = 85;
KSI = Smax^2/Km_fit;

% Part C
y = fit(1)*(1/70) + fit(2);
v = 1/y;
%% Question 4 

% Part b 

% Estimating line fit at drug effect plot 

x = [2 3];
y = [0.5125 0.009168];

plot_fit = polyfit(x,y,1);

% Finding x value for y=0.5 (50% of inhibition) 

IC50_log = (0.5 - plot_fit(2))/plot_fit(1);


