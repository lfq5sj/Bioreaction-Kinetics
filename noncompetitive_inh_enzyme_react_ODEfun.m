
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of ODEs describing enzymatic reaction kinetics in the presence of a competitive inhibitor:
% E + S <-> ES      (k1f, k1r)
% ES -> E + P       (k2f)
% E + I <-> EI      (kif, kir) 
% EI + S <-> ESI    (k1f, k1r)
% ES + I <-> ESI    (kif, kir)


function [dydt] = noncompetitive_inh_enzyme_react_ODEfun(t,y,parameters)

% Assign names for parameters
[k1f, k1r, k2f, kif, kir] = parameters{:};

P = y(1);
ES = y(2);
E = y(3);
I = y(4);
EI = y(5);
S = y(6);
EIS = y(7);

% Differential equations;
dP = k2f*ES;                                                                        % [uM/s] product
dES = k1f*E*S - k1r*ES - k2f*ES -kif*ES*I + kir*EIS;                                % [uM/s] enzyme-substrate complex
dE = -k1f*E*S + k1r*ES + k2f*ES - kif*E*I + kir*EI;                                 % [uM/s] free enzyme
dI = -kif*E*I + kir*EI -kif*ES*I + kir*EIS;                                         % Inhibitor
dEI = kif*E*I - kir*EI - k1f*EI*S + k1r*EIS;                                        % Enzyme-inhibitor complex
dS = -k1f*E*S + k1r*ES - k1f*EI*S + k1r*EIS;                                        % [uM/s] substrate
dEIS = -kir*EIS + kif*ES*I - k1r*EIS + k1f*EI*S;                                    % Enzyme-inhibitor-substrate complex

dydt = [dP;dES;dE;dI;dEI;dS;dEIS];                                                  % Reassemble differential equations