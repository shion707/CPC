% Revised CPC model with DCN -> IO suppression -- natural (non-inverted) variables
% Wang & Ivry (2025), Sci. Adv. 11, eadr4540. Post-hoc variant adding an
% inhibitory pathway from the DCN to the inferior olive (Eqs. 8-9).

%% revised CPC with DCN-IO suppression
%  0 ITI
close all
clear

% --- Model parameters (short-ITI / 0 s; published revised-CPC fit) ---
clear
l       = 1.3;      % PF-PC LTD learning rate l                      (Eq. 2)
f       = 0.5;      % PF-PC forgetting rate f (short ITI)            (Eq. 2)
beta    = 0.1;      % MF-DCN LTP learning rate beta (= paper beta)   (Eq. 3)
alpha   = 0.012;    % MF-DCN forgetting rate alpha                   (Eq. 3)
gamma   = 0.52;     % PC -> DCN scaling factor gamma                 (Eq. 4)
epsilon = 202;      % neural activity -> hand-angle scaling epsilon  (Eq. 5)
omega   = 2.8;      % DCN -> IO suppression strength (short ITI)     (Eqs. 8-9)
w_o     = 1;        % baseline PF-PC synaptic strength w_o           (Eq. 2)
m_o     = 1;        % baseline MF-DCN synaptic strength m_o          (Eq. 3)
m_max   = 1.85;     % maximal MF-DCN synaptic strength m_max         (Eq. 3)
DCN_base = m_o - gamma*w_o;   % baseline DCN activation (no learning), for Eq. 8

rot=[0*zeros(11,1)];
r=[10 10 10 10];
for i=1:18
    rot=[rot; r'];
end
rot=[rot; zeros(50,1)];%
N=length(rot);

kappa=15;                            % von Mises concentration (tuning width)
x=-pi:pi/100:pi-0.000001;
u=-pi:pi/100:pi-0.000001;            % preferred error directions i of the units
tuning=circ_vmpdf(x',u,kappa);       % PC tuning functions (computed, unused below)
v=[cos(u)' sin(u)'];                 % unit tuning-direction vectors v_i  (Eq. 5)
cf_dep=2;                            % (unused)
DCN_list=[];
w=u*0+w_o;           % PF-PC synaptic strength  (starts at w_o, falls via LTD)
m=u*0+m_o;           % MF-DCN synaptic strength (starts at m_o, rises via LTP)
DCN=u*0+DCN_base;    % DCN activation (init at baseline so the Eq. 8 term is defined on trial 1)
for n=1:N-1
    r=rot(n);
    if r>0
        err_dir=0;
        CS=abs(circ_vmpdf(err_dir,u,kappa))./max(circ_vmpdf(0,u,kappa));   % Eq. 1
        % Eqs. 8-9 DCN->IO suppression: cs' = max(0, cs - omega*mean(DCN - DCN_base))
        CS=max(0,CS-omega.*mean(DCN-DCN_base));
    elseif r<0
        err_dir=pi;
        CS=abs(circ_vmpdf(err_dir,u,kappa))./max(circ_vmpdf(0,u,kappa));
        CS=max(0,CS-omega.*mean(DCN-DCN_base));
    else
        CS=0;
    end

    w = w - l.*CS;                                             % PF-PC LTD
    w = w + f.*(w_o-w);                                        % PF-PC LTP
    m = m + beta.*(w_o-w).*(m_max-m) + alpha.*(m_o-m).*(1-CS); % MF-DCN dynamics
    DCN = m - gamma.*w;                                       % DCN activity

    DCN_list(n,:)=DCN;
    m_all(n,:)=m;        % stable MF-DCN synaptic strength, full population (stored, unused)
    vol_all(n,:)=w_o-w;  % volatile PF-PC depression (w_o - w), full population (stored, unused)


    popVec=mean(DCN'.*v);
    [dtheta, drho]=cart2pol(popVec(1),popVec(2));
    hand(n)=drho.*cos(dtheta);                                % Delta hand (Eq. 5)

    mC=mean((m-m_o)'.*v);                                     % stable (DCN) component
    [dtheta, drho]=cart2pol(mC(1),mC(2));
    m_list(n)=drho.*cos(dtheta)*epsilon;
    wC=mean((gamma*(w_o-w))'.*v);                             % volatile (PC) component
    [dtheta, drho]=cart2pol(wC(1),wC(2));
    vol_list(n)=drho.*cos(dtheta)*epsilon;

end

figure
subplot(2,2,1);hold on;
%plot(hand(end-100:end)*epsilon,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(vol_list(10:110),'-k','color',[0.9,0.75,0],'linewidth',2,'markersize',6);    % volatile (PC)
plot(m_list(10:110),'-k','color',[0.3,0.5,0.9],'linewidth',2,'markersize',6);     % stable (DCN)
plot(hand(10:110)*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6); % total hand

axis([-10 110 -5 35])
set(gca, 'LineWidth',1);
set(gca,'xtick',[0:50:800]);
% set(gca,'xticklabel',[0:20:80]*2,'fontsize',15);
set(gca,'ytick',[-20:10:30]);

ax = gca;
ax.YAxis.TickDirection = 'out';
ax.XAxis.TickDirection = 'out';
set(gca, 'LineWidth',1);
box off

hand1=hand;


figure(123)
subplot(2,2,4);hold on;
% plot(-hand(121:N-1)*epsilon,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand(11:110)*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([-10 N+9 -5 35])
box on


% subplot(2,2,2);hold on;
% bfs=[DCN_list(end,:).*circ_vmpdf(x',u,kappa)]';
% for i=1:10:length(u)
%     plot(x,bfs(i,:),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
% end
%
% axis([-pi pi 0 2])
% box on

earlyhand(1)=mean(hand(10:20));
earlywash(1)=power(nanmean(hand(112:121)./nanmean(hand(110))),1);
%%
% 6s ITI

f=0.75;          % PF-PC forgetting rate (long ITI)              (Eq. 2)
omega=0;         % DCN -> IO suppression off in the long-ITI condition (Eqs. 8-9)

x=-pi:pi/100:pi-0.000001;
u=-pi:pi/100:pi-0.000001;
tuning=circ_vmpdf(x',u,kappa);
v=[cos(u)' sin(u)'];
cf_dep=2;
DCN_list=[];
w=u*0+w_o;
m=u*0+m_o;
for n=1:N-1
    r=rot(n);
    if r>0
        err_dir=0;
        CS=abs(circ_vmpdf(err_dir,u,kappa))./max(circ_vmpdf(0,u,kappa));
        CS=max(0,CS-omega.*mean(DCN-DCN_base));
    elseif r<0
        err_dir=pi;
        CS=abs(circ_vmpdf(err_dir,u,kappa))./max(circ_vmpdf(0,u,kappa));
        CS=max(0,CS-omega.*mean(DCN-DCN_base));
    else
        CS=0;
    end

    w = w - l.*CS;
    w = w + f.*(w_o-w);
    m = m + beta.*(w_o-w).*(m_max-m) + alpha.*(m_o-m).*(1-CS);
    DCN = m - gamma.*w;

    DCN_list(n,:)=DCN;
    m_all(n,:)=m;
    vol_all(n,:)=w_o-w;

    popVec=mean(DCN'.*v);
    [dtheta, drho]=cart2pol(popVec(1),popVec(2));
    hand(n)=drho.*cos(dtheta);

    mC=mean((m-m_o)'.*v);
    [dtheta, drho]=cart2pol(mC(1),mC(2));
    m_list(n)=drho.*cos(dtheta)*epsilon;
    wC=mean((gamma*(w_o-w))'.*v);
    [dtheta, drho]=cart2pol(wC(1),wC(2));
    vol_list(n)=drho.*cos(dtheta)*epsilon;

end
%%
figure; hold on

subplot(2,2,1);hold on;
%plot(hand(end-100:end)*epsilon,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(vol_list(10:110),'-k','color',[0.9,0.75,0],'linewidth',2,'markersize',6);    % volatile (PC)
plot(m_list(10:110),'-k','color',[0.3,0.5,0.9],'linewidth',2,'markersize',6);     % stable (DCN)
plot(hand(10:110)*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6); % total hand

axis([-10 110 -5 35])
set(gca, 'LineWidth',1);
set(gca,'xtick',[0:50:800]);
% set(gca,'xticklabel',[0:20:80]*2,'fontsize',15);
set(gca,'ytick',[-20:10:30]);

ax = gca;
ax.YAxis.TickDirection = 'out';
ax.XAxis.TickDirection = 'out';
set(gca, 'LineWidth',1);
box off

%%

figure(123)
%learning curve
subplot(2,2,4);hold on;
% plot(-hand(121:N-1)*epsilon,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand(11:110)*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([-10 110 -5 35])
box on

%ratio
subplot(2,2,3);hold on;
% plot(-hand(121:N-1)*epsilon,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand1(11:110)./hand(11:110),'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([-10 110 0 4])
box on