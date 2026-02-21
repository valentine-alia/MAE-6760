%% MAE 6760 Model Based Estimation
% Cornell University
% M Campbell
%
% Homework #2
% Problem #4: MMSE Estimation
%   for a time varying output function
%
clear all;close all

%% simulate the measurements
t=[0.1:0.1:1]';
n=length(t);
x=[1;1;1];
R=diag([0.001 0.002  0.005  0.010  0.008  0.002  0.010  0.007  0.020 0.006].^2);
v=sqrtm(R)*randn(n,1);
z=x(1)+x(2)*sin(10*t)+x(3)*exp(2*t.^2)+v;
z

%% prior information
x0=[1.01;0.98;0.99];
P0=eye(3)*0.001;

%below is a plotting script with an error bar if it helps...
%alphai = an (ni x 1) vector of multiplicative factors
%xi = 3 state estimates as a function of alphai (ni x 3)
%sigi = 3 standard deviationsn for each state estimates as a function of alphai (ni x 3)
your_title='fill this in...';
figure;
t = tiledlayout(1,3,'TileSpacing','compact'); % 1 row, 3 columns
nexttile;
errorbar(alphai,xi(1,:),sigi(1,:)*2,'b');
set(gca,'XScale','log','xtick',10.^[-3 -2 -1 0 1 2 3]);grid;
hold on;semilogx(alphai,ones(ni,1)*1,'color',[200 0 0]/255);hold off;
xlabel('scale factor \alpha');ylabel('state estimate');
nexttile;
errorbar(alphai,xi(2,:),sigi(2,:)*2,'b');
set(gca,'XScale','log','xtick',10.^[-3 -2 -1 0 1 2 3]);grid;
hold on;semilogx(alphai,ones(ni,1)*1,'color',[200 0 0]/255);hold off;
legend('est +/- 2\sigma','true value','Location','south');
ylabel('state estimate');
xlabel('scale factor \alpha');
nexttile;
errorbar(alphai,xi(3,:),sigi(3,:)*2,'b');
set(gca,'XScale','log','xtick',10.^[-3 -2 -1 0 1 2 3]);grid;
hold on;semilogx(alphai,ones(ni,1)*1,'color',[200 0 0]/255);hold off;
xlabel('scale factor \alpha');ylabel('state estimate');
title(t,your_title);
PrepTilePresentation(gcf);%makes lines thicker and font bigger


%% Extra functions for plotting
function PrepTilePresentation(fig_num);
%
% This function prepares a figure for presentations
% with thicker lines and bigger font
%   Fontsize: 12
%   Fontweight: bold
%   LineWidth: 2
% 
% INPUTS
%  fig_num = figure number 
%
fig_handle=figure(fig_num);
fig_children=get(fig_handle,'children'); %find all sub-plots
fig_children_children=get(fig_children,'Children');

set(fig_children.Title,'FontSize',12,'FontWeight','bold');
set(fig_children.XLabel,'FontSize',12,'FontWeight','bold');
set(fig_children.YLabel,'FontSize',12,'FontWeight','bold');

for i=1:length(fig_children_children),
    
    set(fig_children_children(i),'FontSize',12);
    set(fig_children_children(i),'FontWeight','bold');
    
    fig_children_children_children=get(fig_children_children(i),'Children');
    set(fig_children_children_children,'LineWidth',2);
end

end