%% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
% Homework #2
% Problem #2: Parameter Estimation (continuous)
%   using a triangular non-Gaussian noise distribution
%

%% define triangular noise pdf
pdv = makedist('Triangular','a',0,'b',0,'c',1.0);
mn_v=mean(pdv); %mean of the noise distribution
sig_v=std(pdv); %std of the noise distribution
disp('mean and std of the sensor noise are:')
disp([mn_v sig_v]);
v = [-1.5:0.001:1.5];
fv = pdf(pdv,v);
%plot pdf
figure(1);
plot(v,fv,'b-');
ylabel('Sensor Noise Model f(v)');
xlabel('v');
grid;
PrepFigPresentation(1);

%% define uniform prior pdfsize(
pdx0 = makedist('Uniform','Lower',3,'Upper',4);
mn_x0=mean(pdx0); %mean of the prior distribution
sig_x0=std(pdx0); %std of the noise distribution
disp('mean and std of the prior are:')
disp([mn_x0 sig_x0]);
x = [2:0.001:5];
fx0 = pdf(pdx0,x);
%plot pdf
figure(2);
plot(x,fx0,'b-');
ylabel('Prior Distribution f(x_0)');
xlabel('x_0');
grid;
PrepFigPresentation(2);


%% Extra functions for plotting
function PrepFigPresentation(fignum);
%
% prepares a figure for presentations
%
% Fontsize: 16
% Fontweight: bold
% LineWidth: 2
% 
figure(fignum);
fig_children=get(fignum,'children'); %find all sub-plots

for i=1:length(fig_children),
    
    set(fig_children(i),'FontSize',16);
    set(fig_children(i),'FontWeight','bold');
    
    fig_children_children=get(fig_children(i),'Children');
    set(fig_children_children,'LineWidth',2);
end
end