% CPC model
% Wang & Ivry (2025), Sci. Adv. 11, eadr4540.

%   Eq. 1  CS_i = VM(theta_e, i, s)
%   Eq. 2  w update (LTD + forgetting toward w_o)
%   Eq. 3  m update (LTP gated by the PF-PC depression w_o - w, bounded by m_max)
%   Eq. 4  DCN_i = m - gamma*w
%   Eq. 5  hand  = population resultant of v_i*DCN_i (projected on the error axis)


%% Exp 1
close all
clear

% --- Model parameters (fitted to Exp 1; held fixed across all experiments) ---
l       = 2;        % PF-PC LTD learning rate l                      (Eq. 2)
f       = 0.5;      % PF-PC forgetting rate f                        (Eq. 2)
beta    = 0.05;     % MF-DCN LTP learning rate beta                  (Eq. 3)
alpha   = 0.018;    % MF-DCN forgetting rate alpha                   (Eq. 3)
gamma   = 0.15;     % PC -> DCN scaling factor gamma                 (Eq. 4)
epsilon = 130;      % neural activity -> hand-angle scaling epsilon  (Eq. 5)
w_o     = 1;        % baseline PF-PC synaptic strength w_o           (Eq. 2)
m_o     = 1;        % baseline MF-DCN synaptic strength m_o          (Eq. 3)
m_max   = 1.85;     % maximal MF-DCN synaptic strength m_max         (Eq. 3)

% --- Perturbation (clamp) schedule, one entry per trial ---
rot = [0*ones(30,1); 10*ones(100,1); 0*ones(70,1)];
N   = length(rot);

% --- Population of tuned units ---
kappa  = 15;                        % von Mises concentration (tuning width)
x      = -pi:pi/1000:pi;
u      = -pi:pi/1000:pi;            % preferred error directions i of the units
tuning = circ_vmpdf(x',u,kappa);    % PC tuning functions (computed, unused below)
v      = [cos(u)' sin(u)'];         % unit tuning-direction vectors v_i      (Eq. 5)

% --- Learning ---
DCN_list = [];
w = u*0 + w_o;       % PF-PC synaptic strength  (starts at w_o, falls via LTD)
m = u*0 + m_o;       % MF-DCN synaptic strength (starts at m_o, rises via LTP)
for n = 1:N-1
    r = rot(n);                                  % perturbation on trial n
    if r > 0
        err_dir = 0;                             % error direction theta_e
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));   % Eq. 1
    elseif r < 0
        err_dir = pi;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    else
        CS = 0;
    end

    % PF-PC LTD: complex spikes depress the synapse
    w = w - l.*CS;

    % MF-DCN LTP gated by the current PF-PC depression (w_o - w),
    %         bounded by m_max, plus passive forgetting toward m_o
    m = m + beta.*(w_o - w).*(m_max - m) + alpha.*(m_o - m);

    % PF-PC passive forgetting toward w_o
    w = w + f.*(w_o - w);

    % DCN activation per unit
    DCN = m - gamma.*w;

    % population read-out: resultant of v_i*DCN_i, projected on the error
    %         axis (positive = adaptation, opposite the clamp)
    DCN_list(n,:) = DCN;
    popVec = mean(DCN'.*v);
    [dtheta, drho] = cart2pol(popVec(1),popVec(2));
    hand(n) = drho.*cos(dtheta);                 % Delta hand (Eq. 5, before epsilon)

    % Decomposition of the hand angle into its two contributions (for plotting)
    mC = mean((m-m_o)'.*v);                      % stable (DCN/m) component
    [dtheta, drho] = cart2pol(mC(1),mC(2));
    m_list(n) = drho.*cos(dtheta)*epsilon;
    wC = mean((gamma*(w_o-w))'.*v);              % volatile (PC) component, gamma*(w_o-w)
    [dtheta, drho] = cart2pol(wC(1),wC(2));
    vol_list(n) = drho.*cos(dtheta)*epsilon;
end

figure
subplot(2,2,1);hold on;
plot(vol_list,'-k','color',[0.9,0.75,0],'linewidth',2,'markersize',6);      % volatile (PC / cerebellar cortex)
plot(m_list,'-k','color',[0.3,0.5,0.9],'linewidth',2,'markersize',6);       % stable (DCN)
plot(hand*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6); % total hand angle
axis([0 N -5 25])


%% Exp 2
close all
clear

% --- Model parameters (same fixed set as Exp 1; Eqs. 1-5) ---
l       = 2;        % PF-PC LTD learning rate l                      (Eq. 2)
f       = 0.5;      % PF-PC forgetting rate f                        (Eq. 2)
beta    = 0.05;     % MF-DCN LTP learning rate beta                  (Eq. 3)
alpha   = 0.018;    % MF-DCN forgetting rate alpha                   (Eq. 3)
gamma   = 0.15;     % PC -> DCN scaling factor gamma                 (Eq. 4)
epsilon = 130;      % neural activity -> hand-angle scaling epsilon  (Eq. 5)
w_o     = 1;        % baseline PF-PC synaptic strength w_o           (Eq. 2)
m_o     = 1;        % baseline MF-DCN synaptic strength m_o          (Eq. 3)
m_max   = 1.85;     % maximal MF-DCN synaptic strength m_max         (Eq. 3)

% --- Perturbation (clamp) schedule ---
rot = [0*ones(10,1); 10*ones(100,1)];
for i = 1:30
    rot = [rot; 0; 10];
end
N = length(rot);

% --- Population of tuned units ---
kappa  = 15;                        % von Mises concentration (tuning width)
x      = -pi:pi/1000:pi;
u      = -pi:pi/1000:pi;            % preferred error directions i of the units
tuning = circ_vmpdf(x',u,kappa);    % PC tuning functions (computed, unused below)
v      = [cos(u)' sin(u)'];         % unit tuning-direction vectors v_i      (Eq. 5)

% --- Learning ---
DCN_list = [];
w = u*0 + w_o;       % PF-PC synaptic strength  (starts at w_o, falls via LTD)
m = u*0 + m_o;       % MF-DCN synaptic strength (starts at m_o, rises via LTP)
for n = 1:N-1
    r = rot(n);
    if r > 0
        err_dir = 0;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    elseif r < 0
        err_dir = pi;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    else
        CS = 0;
    end

    w = w - l.*CS;                                         % Eq. 2a LTD
    m = m + beta.*(w_o-w).*(m_max-m) + alpha.*(m_o-m);     % Eq. 3 LTP + forgetting
    w = w + f.*(w_o-w);                                    % Eq. 2b forgetting
    DCN = m - gamma.*w;                                    % Eq. 4
%     m_list(n,:)=m;
%     vol_list(n,:)=w_o-w;
    DCN_list(n,:) = DCN;
    popVec = mean(DCN'.*v);
    [dtheta, drho] = cart2pol(popVec(1),popVec(2));
    hand(n) = drho.*cos(dtheta);                           % Delta hand (Eq. 5)

    mC = mean((m-m_o)'.*v);                                % stable (DCN/m) component
    [dtheta, drho] = cart2pol(mC(1),mC(2));
    m_list(n) = drho.*cos(dtheta)*epsilon;
    wC = mean((gamma*(w_o-w))'.*v);                        % volatile (PC) component
    [dtheta, drho] = cart2pol(wC(1),wC(2));
    vol_list(n) = drho.*cos(dtheta)*epsilon;
end

% figure; hold on
% subplot(3,2,1);hold on;
% bs_circle=DCN_list(end,:)'.*v;
% % plot(bs_circle(:,2),bs_circle(:,1),'-k','color',[0.5,0.8,0.7],'linewidth',2);
% fill(sin(x),cos(x),[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.1);
% for i=1:length(u)-1
%     plot(bs_circle(i:i+1,2),bs_circle(i:i+1,1),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',3,'markersize',20);
% end
% axis([-2 2 -2 2],'square')
%
% subplot(3,2,3);hold on;
% plot(hand*epsilon,'-k');
%

figure
subplot(2,2,1);hold on;
plot(vol_list,'-k','color',[0.9,0.75,0],'linewidth',2,'markersize',6);      % volatile (PC / cerebellar cortex)
plot(m_list,'-k','color',[0.3,0.5,0.9],'linewidth',2,'markersize',6);       % stable (DCN)
plot(hand*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6); % total hand angle
axis([0 N -5 25])

%% Exp 3 -- Anterograde interference
close all
clear

% --- Model parameters (same fixed set as Exp 1; Eqs. 1-5) ---
l       = 2;        % PF-PC LTD learning rate l                      (Eq. 2)
f       = 0.5;      % PF-PC forgetting rate f                        (Eq. 2)
beta    = 0.05;     % MF-DCN LTP learning rate beta                  (Eq. 3)
alpha   = 0.018;    % MF-DCN forgetting rate alpha                   (Eq. 3)
gamma   = 0.15;     % PC -> DCN scaling factor gamma                 (Eq. 4)
epsilon = 130;      % neural activity -> hand-angle scaling epsilon  (Eq. 5)
w_o     = 1;        % baseline PF-PC synaptic strength w_o           (Eq. 2)
m_o     = 1;        % baseline MF-DCN synaptic strength m_o          (Eq. 3)
m_max   = 1.85;     % maximal MF-DCN synaptic strength m_max         (Eq. 3)

% --- Perturbation (clamp) schedule: 30 clamp reversed for 200 trials ---
rot = [0*ones(10,1); 10*ones(100,1); -10*ones(200,1)];
N = length(rot);

% --- Population of tuned units ---
kappa  = 15;                        % von Mises concentration (tuning width)
x      = -pi:pi/100:pi;
u      = -pi:pi/100:pi;             % preferred error directions i of the units
tuning = circ_vmpdf(x',u,kappa);    % PC tuning functions (computed, unused below)
v      = [cos(u)' sin(u)'];         % unit tuning-direction vectors v_i      (Eq. 5)

% --- Learning ---
DCN_list = [];
w = u*0 + w_o;       % PF-PC synaptic strength  (starts at w_o, falls via LTD)
m = u*0 + m_o;       % MF-DCN synaptic strength (starts at m_o, rises via LTP)
for n = 1:N-1
    r = rot(n);
    if r > 0
        err_dir = 0;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    elseif r < 0
        err_dir = pi;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    else
        CS = 0;
    end

    w = w - l.*CS;                                         % Eq. 2a LTD
    m = m + beta.*(w_o-w).*(m_max-m) + alpha.*(m_o-m);     % Eq. 3 LTP + forgetting
    w = w + f.*(w_o-w);                                    % Eq. 2b forgetting
    DCN = m - gamma.*w;                                    % Eq. 4
    m_list(n,:) = m;            % stable MF-DCN synaptic strength, full population
    vol_list(n,:) = w_o - w;    % volatile PF-PC depression (w_o - w), full population
    DCN_list(n,:) = DCN;
    popVec = mean(DCN'.*v);
    [dtheta, drho] = cart2pol(popVec(1),popVec(2));
    hand(n) = drho.*cos(dtheta);                           % Delta hand (Eq. 5)
end

figure; hold on

figure
subplot(2,2,1);hold on;

plot(hand*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6); % total hand angle
axis([0 N -25 25])

%% Exp 5 -- p(switch) = 0.125
close all
clear

% --- Model parameters (same fixed set as Exp 1; Eqs. 1-5) ---
l       = 2;        % PF-PC LTD learning rate l                      (Eq. 2)
f       = 0.5;      % PF-PC forgetting rate f                        (Eq. 2)
beta    = 0.05;     % MF-DCN LTP learning rate beta                  (Eq. 3)
alpha   = 0.018;    % MF-DCN forgetting rate alpha                   (Eq. 3)
gamma   = 0.15;     % PC -> DCN scaling factor gamma                 (Eq. 4)
epsilon = 130;      % neural activity -> hand-angle scaling epsilon  (Eq. 5)
w_o     = 1;        % baseline PF-PC synaptic strength w_o           (Eq. 2)
m_o     = 1;        % baseline MF-DCN synaptic strength m_o          (Eq. 3)
m_max   = 1.85;     % maximal MF-DCN synaptic strength m_max         (Eq. 3)

tim = 1:10;
yt = 0*[tim tim tim tim tim tim tim tim];

rot = [0*ones(10,1); 10*ones(100,1); -10*ones(12,1); 0*ones(30,1)];

% change here for the switch probability
r = [10 10 10 10 10 10 10 10 -10 -10 -10 -10 -10 -10 -10 -10 ]; % p(switch) = 0.125

for i = 1:20
    rot = [rot; r'];
end

rot = [rot; 0*ones(10,1); 10*ones(200,1)];

N = length(rot);

% --- Population of tuned units ---
kappa  = 15;                        % von Mises concentration (tuning width)
x      = -pi:pi/100:pi;
u      = -pi:pi/100:pi;             % preferred error directions i of the units
tuning = circ_vmpdf(x',u,kappa);    % PC tuning functions (computed, unused below)
v      = [cos(u)' sin(u)'];         % unit tuning-direction vectors v_i      (Eq. 5)

% plot(x,y)
figure; hold on
set(0,'defaultfigurecolor','w')
% --- Learning ---
DCN_list = [];
w = u*0 + w_o;       % PF-PC synaptic strength  (starts at w_o, falls via LTD)
m = u*0 + m_o;       % MF-DCN synaptic strength (starts at m_o, rises via LTP)
for n = 1:N-1
    r = rot(n);
    if r > 0
        err_dir = 0;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    elseif r < 0
        err_dir = pi;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    else
        CS = 0;
    end

    w = w - l.*CS;                                         % Eq. 2a LTD
    m = m + beta.*(w_o-w).*(m_max-m) + alpha.*(m_o-m);     % Eq. 3 LTP + forgetting
    w = w + f.*(w_o-w);                                    % Eq. 2b forgetting
    DCN = m - gamma.*w;                                    % Eq. 4
    m_list(n,:) = m;            % stable MF-DCN synaptic strength, full population
    vol_list(n,:) = w_o - w;    % volatile PF-PC depression (w_o - w), full population
    DCN_list(n,:) = DCN;
    popVec = mean(DCN'.*v);
    [dtheta, drho] = cart2pol(popVec(1),popVec(2));
    hand(n) = drho.*cos(dtheta);                           % Delta hand (Eq. 5)
end

figure; hold on

% %PC in polar
% subplot(2,2,1);hold on;
% bs_circle=DCN_list(end,:)'.*v;
% % plot(bs_circle(:,2),bs_circle(:,1),'-k','color',[0.5,0.8,0.7],'linewidth',2);
% fill(sin(x),cos(x),[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.1);
% for i=1:length(u)-1
%     plot(bs_circle(i:i+1,2),bs_circle(i:i+1,1),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',3,'markersize',20);
% end
% axis([-2 2 -2 2],'square')
%
% %PC
% subplot(2,2,2);hold on;
% bfs=[DCN_list(end,:).*circ_vmpdf(x',u,kappa)]';
% for i=1:10:length(u)
%     plot(x,bfs(i,:),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
% end
%
% axis([-pi pi 0 2])
%
%
%
% % learning curve
% subplot(2,2,3);hold on;
% plot(hand*epsilon,'.k');
% axis([0 N -30 30])
% learning vs relearning
subplot(2,2,4);hold on;
plot(hand(end-200:end)*epsilon,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand(10:110)*epsilon,'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([0 100 -5 30])

%% TBT (trial-by-trial) design
close all
clear

% --- Model parameters (same fixed set as Exp 1; Eqs. 1-5) ---
l       = 2;        % PF-PC LTD learning rate l                      (Eq. 2)
f       = 0.5;      % PF-PC forgetting rate f                        (Eq. 2)
beta    = 0.05;     % MF-DCN LTP learning rate beta                  (Eq. 3)
alpha   = 0.018;    % MF-DCN forgetting rate alpha                   (Eq. 3)
gamma   = 0.15;     % PC -> DCN scaling factor gamma                 (Eq. 4)
epsilon = 130;      % neural activity -> hand-angle scaling epsilon  (Eq. 5)
w_o     = 1;        % baseline PF-PC synaptic strength w_o           (Eq. 2)
m_o     = 1;        % baseline MF-DCN synaptic strength m_o          (Eq. 3)
m_max   = 1.85;     % maximal MF-DCN synaptic strength m_max         (Eq. 3)

tim = 1:10;

yt = 0*[tim tim tim tim tim tim tim tim];

rot=[]

r = [10 -10 10 -10 ];
for i = 1:300
    rot = [rot; r(randperm(4))'];
end


N = length(rot);

% --- Population of tuned units ---
kappa  = 15;                        % von Mises concentration (tuning width)
x      = -pi:pi/100:pi;
u      = -pi:pi/100:pi;             % preferred error directions i of the units
tuning = circ_vmpdf(x',u,kappa);    % PC tuning functions (computed, unused below)
v      = [cos(u)' sin(u)'];         % unit tuning-direction vectors v_i      (Eq. 5)

% --- Learning ---
DCN_list = [];
w = u*0 + w_o;       % PF-PC synaptic strength  (starts at w_o, falls via LTD)
m = u*0 + m_o;       % MF-DCN synaptic strength (starts at m_o, rises via LTP)
for n = 1:N-1
    r = rot(n);
    if r > 0
        err_dir = 0;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    elseif r < 0
        err_dir = pi;
        CS = circ_vmpdf(err_dir,u,kappa)./max(circ_vmpdf(0,u,kappa));
    else
        CS = 0;
    end

    w = w - l.*CS;                                         % Eq. 2a LTD
    m = m + beta.*(w_o-w).*(m_max-m) + alpha.*(m_o-m);     % Eq. 3 LTP + forgetting
    w = w + f.*(w_o-w);                                    % Eq. 2b forgetting

    DCN = m - gamma.*w;                                    % Eq. 4

    m_list(n,:) = m;            % stable MF-DCN synaptic strength, full population
    vol_list(n,:) = w_o - w;    % volatile PF-PC depression (w_o - w), full population
    DCN_list(n,:) = DCN;

    popVec = mean(DCN'.*v);
    [dtheta, drho] = cart2pol(popVec(1),popVec(2));
    hand(n) = drho.*cos(dtheta);                           % Delta hand (Eq. 5)
end


figure(114)
k=0;
for n=[1 2 3 4 50 51 53 54]
    k=k+1;
    if k<5
        subplot(2,2,1);hold on;
    else
        subplot(2,2,2);hold on;
    end
    bfs=[m_list(n,:)-m_o]';        % stable LTP gain (m - m_o), full population
    plot(x,bfs,'-r','color',[0.7,0.8-0.2*mod(n,50),1-0.1*mod(n,50)],'linewidth',1.5);

    % for i=1:length(u)
    %     plot(x(i),bfs(i),'.r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
    % end
    box off
    axis([-pi pi 0 1])
    set(gca, 'LineWidth',1);
    % axis([0.5 3.5 -2 4])
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
end


% Fig S5
figure(115)
k=0;
for n=[1 2 3 4 51 52 53 54]
    k=k+1;
    if k<5
        subplot(2,2,1);hold on;
    else
        subplot(2,2,2);hold on;
    end
    bfs=[2-vol_list(n,:)]';        % vol_list = depression (w_o - w)
    plot(x,bfs,'-r','color',[0.7,0.8-0.2*mod(n,50),1-0.1*mod(n,50)],'linewidth',1.5);

    % for i=1:length(u)
    %     plot(x(i),bfs(i),'.r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
    % end
    box off
    axis([-pi pi 0 4])
    set(gca, 'LineWidth',1);
    % axis([0.5 3.5 -2 4])
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
end


dccw=find(rot>0);
dcw=find(rot<0);
hand=hand*epsilon;
hand(end+1:end+5)=nan;
dhand=hand(2:end)-hand(1:end-1);
noflip=sign(rot(1:end-1).*rot(2:end));

dhand2=dhand;
dhand2(dccw)=-dhand2(dccw);
dhand2=dhand2(200:1000);
nanmean(dhand2)

dhand=[0,dhand];
dhand1=dhand;
dhand1(dccw)=-dhand1(dccw);
% dhand1(1)=[];
% dhnof=dhand1((noflip>0));
% dhf=dhand1((noflip<0));
dhandTS=nanmean(reshape(dhand1(1:200),20,10));

dhand1=dhand1(200:1000);
nanmean(dhand1)
% nanmean(dhf)
% nanmean(dhnof)

dhand3=dhand(3:end);
dhand3(dccw)=-dhand3(dccw);
dhand3=dhand3(200:1000);
nanmean(dhand3);

%%
% Exp 5 50% condition
figure
subplot(2,3,1);hold on;
plot([-100 10000],[0 0],'k--','linewidth',0.5)
h=-[dhand1;dhand2;dhand3];
plot(nanmean(h'),'.-','color',[0.2 0.5 1],'linewidth',2,'markersize',20)
set(gca,'xtick',[]);
axis([0 4 -5 5])


%retention
subplot(2,3,2);hold on;
h=1+nanmean(dhand2)./nanmean(dhand1);
bar(nanmean(h'))
set(gca,'xtick',[]);
axis([0 2 0 1])

%learning
subplot(2,3,3);hold on;
h=-[dhand1];
bar(nanmean(h'))
set(gca,'xtick',[]);
axis([0 4 0 5])


% subplot(2,3,4);hold on;
% plot(-dhandTS,'-','color',[0.2 0.5 1],'linewidth',2,'markersize',20)
% set(gca,'xtick',[]);
% %axis([0 10 0 5])