%% revised CPC with DCN-IO surpression
% Exp 10: 0 ITI
close all
clear

%initial pramater

clear
w0=1;
w_lowb=0.15;
r_ltd=0.2;
r_ltp=0.015;
fCS=0.5;
PCd=0.5;
PCa=2;
scale1=0.2;
scale2=210;
scale3=2.5;

rot=[0*zeros(11,1)];
r=[10 10 10 10];
for i=1:18
    rot=[rot; r'];
end
rot=[rot; zeros(50,1)];%
N=length(rot);

sd=15;
x=-pi:pi/100:pi-0.000001;
u=-pi:pi/100:pi-0.000001;
PC=circ_vmpdf(x',u,sd);
basevec=[cos(u)' sin(u)'];
cf_dep=2;
ss_list=[];
w=u*0+w0;
cs_PC=u*0;
ss_PC=0;
for n=1:N-1
    r=rot(n);
    if r>0
        error=0;
        act_CF=abs(circ_vmpdf(error,u,sd))./max(circ_vmpdf(0,u,sd));
        act_CF=max(0,act_CF-scale3.*mean(1-ss_PC));
    elseif r<0
        error=pi;
        act_CF=abs(circ_vmpdf(error,u,sd))./max(circ_vmpdf(0,u,sd));
        act_CF=max(0,act_CF-scale3.*mean(1-ss_PC));
    else
        act_CF=0;
    end
    
    cs_PC=cs_PC+(PCa).*act_CF;
    cs_PC=cs_PC-PCd.*cs_PC;
    w=w-r_ltp.*(w-w0).*(1-(act_CF))-cs_PC.*r_ltd.*fCS.*(w-w_lowb);
    ss_PC=w-cs_PC*scale1;
    
    ss_list(n,:)=ss_PC;
    wl(n,:)=w;
    csl(n,:)=cs_PC;
    
    
    DCN=mean(ss_PC'.*basevec);
    [dtheta, drho]=cart2pol(DCN(1),DCN(2));
    hand(n)=-drho.*cos(dtheta);
    
    wT=mean(w'.*basevec);
    [dtheta, drho]=cart2pol(wT(1),wT(2));
    w_list(n)=-drho.*cos(dtheta)*scale2;
    csT=mean(cs_PC'.*basevec.*scale1);
    [dtheta, drho]=cart2pol(csT(1),csT(2));
    cs_list(n)=-drho.*cos(dtheta)*scale2;
    
end

figure
subplot(2,2,1);hold on;
%plot(hand(end-100:end)*scale2,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(-cs_list(10:110),'-k','color',[0.9,0.75,0],'linewidth',2,'markersize',6);
plot(w_list(10:110),'-k','color',[0.3,0.5,0.9],'linewidth',2,'markersize',6);
plot(hand(10:110)*scale2,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6);

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
% plot(-hand(121:N-1)*scale2,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand(11:110)*scale2,'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([-10 N+9 -5 35])
box on


% subplot(2,2,2);hold on;
% bfs=[ss_list(end,:).*circ_vmpdf(x',u,sd)]';
% for i=1:10:length(u)
%     plot(x,bfs(i,:),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
% end
% 
% axis([-pi pi 0 2])
% box on

earlyhand(1)=mean(hand(10:20));
earlywash(1)=power(nanmean(hand(112:121)./nanmean(hand(110))),1);
%%
% Exp 10: 6s ITI

PCd=0.75;
scale3=0;

x=-pi:pi/100:pi-0.000001;
u=-pi:pi/100:pi-0.000001;
PC=circ_vmpdf(x',u,sd);
basevec=[cos(u)' sin(u)'];
cf_dep=2;
ss_list=[];
w=u*0+w0;
cs_PC=u*0;
for n=1:N-1
    r=rot(n);
    if r>0
        error=0;
        act_CF=abs(circ_vmpdf(error,u,sd))./max(circ_vmpdf(0,u,sd));
        act_CF=max(0,act_CF-scale3.*mean(1-ss_PC));
    elseif r<0
        error=pi;
        act_CF=abs(circ_vmpdf(error,u,sd))./max(circ_vmpdf(0,u,sd));
        act_CF=max(0,act_CF-scale3.*mean(1-ss_PC));
    else
        act_CF=0;
    end
    
    cs_PC=cs_PC+(PCa).*act_CF;
    cs_PC=cs_PC-PCd.*cs_PC;
    w=w-r_ltp.*(w-w0).*(1-(act_CF))-cs_PC.*r_ltd.*fCS.*(w-w_lowb);
    ss_PC=w-cs_PC*scale1;
    
    ss_list(n,:)=ss_PC;
    wl(n,:)=w;
    csl(n,:)=cs_PC;
    
    DCN=mean(ss_PC'.*basevec);
    [dtheta, drho]=cart2pol(DCN(1),DCN(2));
    hand(n)=-drho.*cos(dtheta);
    
    wT=mean(w'.*basevec);
    [dtheta, drho]=cart2pol(wT(1),wT(2));
    w_list(n)=-drho.*cos(dtheta)*scale2;
    csT=mean(cs_PC'.*basevec.*scale1);
    [dtheta, drho]=cart2pol(csT(1),csT(2));
    cs_list(n)=-drho.*cos(dtheta)*scale2;
    
end
%%
figure; hold on

subplot(2,2,1);hold on;
%plot(hand(end-100:end)*scale2,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(-cs_list(10:110),'-k','color',[0.9,0.75,0],'linewidth',2,'markersize',6);
plot(w_list(10:110),'-k','color',[0.3,0.5,0.9],'linewidth',2,'markersize',6);
plot(hand(10:110)*scale2,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6);

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
% plot(-hand(121:N-1)*scale2,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand(11:110)*scale2,'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([-10 110 -5 35])
box on

%ratio
subplot(2,2,3);hold on;
% plot(-hand(121:N-1)*scale2,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand1(11:110)./hand(11:110),'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([-10 110 0 4])
box on






%%
% subplot(2,2,2);hold on;
% bfs=[ss_list(end,:).*circ_vmpdf(x',u,sd)]';
% for i=1:10:length(u)
%     plot(x,bfs(i,:),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
% end
% 
% axis([-pi pi 0 2])
% box on

% earlyhand(2)=mean(hand(11:21));
% earlywash(2)=power(nanmean(hand(111:120)./nanmean(hand(90:110))),1);


% figure; subplot(2,2,1)
% %early learning
% plot(earlyhand*scale2,'o-k','color',[0.4,0.6,0.6],'linewidth',2);
% axis([0 3 -5 25])
% set(gca, 'LineWidth',1);
% set(gca,'xtick',[]);
% % set(gca,'xticklabel',[0:20:80]*2,'fontsize',15);
% set(gca,'ytick',[-20:10:30]);
% 
% ax = gca;
% ax.YAxis.TickDirection = 'out';
% ax.XAxis.TickDirection = 'out';
% set(gca, 'LineWidth',1);
% box off
% 
% %early washout
% subplot(2,2,2)
% plot(earlywash,'o-k','color',[0.4,0.6,0.6],'linewidth',2);
% axis([0 3 0.5 1.])
% set(gca, 'LineWidth',1);
% set(gca,'xtick',[]);
% % set(gca,'xticklabel',[0:20:80]*2,'fontsize',15);
% set(gca,'ytick',[-20:0.05:30]);
% 
% ax = gca;
% ax.YAxis.TickDirection = 'out';
% ax.XAxis.TickDirection = 'out';
% set(gca, 'LineWidth',1);
% box off