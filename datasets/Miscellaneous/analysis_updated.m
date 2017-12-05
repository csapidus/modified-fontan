%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Mahdi Al-Husseini, Dataset Analysis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data Import

inlet_time = csvread('inlet_original.csv');
inlet_time = inlet_time(:,1);
inlet_pressure = csvread('inlet_original.csv',0,1,[0,1,137-1,1]) * 133.322365;

LPAOutlet_time = csvread('LPAoutlet_original.csv');
LPAOutlet_time = LPAOutlet_time(:,1);
LPAOutlet_pressure = csvread('LPAoutlet_original.csv',0,1,[0,1,118-1,1]) * 133.322365;

RPAOutlet_time = csvread('RPAoutlet_original.csv');
RPAOutlet_time = RPAOutlet_time(:,1);
RPAOutlet_pressure = csvread('RPAoutlet_original.csv',0,1,[0,1,118-1,1]) * 133.322365;

%% Data Analysis

x = 0.0:0.0001:0.5365;

inlet = pchip(inlet_time, inlet_pressure,x);
LPA = pchip(LPAOutlet_time, LPAOutlet_pressure,x);
RPA = pchip(RPAOutlet_time, RPAOutlet_pressure,x);




%% Data Plotting
plot(x,inlet,'b')
hold on
plot(x,LPA,'r')
hold on
plot(x,RPA,'g')
hold on
legend('Inlet', 'LPA', 'RPA')
xlabel('time (s)')
ylabel('pressure (kPA)')

newInlet(:, 1) = x; 
newInlet(:, 2) = inlet; 
newLPA (:, 1) = x;
newLPA(:, 2) = LPA; 
newRPA (:, 1) = x;
newRPA(:, 2) = RPA; 


%% Saving interpolated data to new CSV

csvwrite('Inlet.csv', newInlet)
csvwrite('Left.csv',newLPA)
csvwrite('Right.csv',newRPA)