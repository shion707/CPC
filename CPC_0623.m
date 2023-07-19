% code of the CPC model for all Exps
%%
%Exp 1

close all
clear

%initial pramater
error_d=0;
w0=1;
w_lowb=0.15;
r_ltd=0.1;
r_ltp=0.018;
fCS=0.5;
PCd=0.5;
PCa=2;
scale1=0.15;
scale2=130;
% pertubation schedule
rot=[0*ones(10,1); 10*ones(100,1);0*ones(60,1)];
N=length(rot);

%base functions
sd=15;
x=-pi:pi/100:pi;
u=-pi:pi/100:pi;
PC=circ_vmpdf(x',u,sd);
basevec=[cos(u)' sin(u)'];

%learning
ss_list=[];
w=u*0+w0;
cs_PC=u*0;
for n=1:N-1
    r=rot(n);
    if r>0
        error=0;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    elseif r<0
        error=pi;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    else
        act_CF=0;
    end
    
    cs_PC=cs_PC+PCa.*act_CF;
    
    w=w-r_ltp.*(w-w0)-cs_PC.*r_ltd.*(w-w_lowb).*fCS;
    cs_PC=cs_PC-PCd.*cs_PC;
    ss_PC=w-cs_PC*scale1;
    w_list(n,:)=w;
    cs_list(n,:)=cs_PC;
    ss_list(n,:)=ss_PC;
    DCN=mean(ss_PC'.*basevec);
    [dtheta, drho]=cart2pol(DCN(1),DCN(2));
    hand(n)=-drho.*cos(dtheta);
    
    
end

figure; hold on
% subplot(3,2,1);hold on;
% bs_circle=ss_list(end,:)'.*basevec;
% % plot(bs_circle(:,2),bs_circle(:,1),'-k','color',[0.5,0.8,0.7],'linewidth',2);
% fill(sin(x),cos(x),[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.1);
% for i=1:length(u)-1
%     plot(bs_circle(i:i+1,2),bs_circle(i:i+1,1),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',3,'markersize',20);
% end
% axis([-2 2 -2 2],'square')

subplot(3,2,1);hold on;
plot(hand*scale2,'-k');
axis([0 N -5 25])

%%
%Exp 2
%initial pramater

close all
clear

%initial pramater
error_d=0;
w0=1;
w_lowb=0.15;
r_ltd=0.1;
r_ltp=0.018;
fCS=0.5;
PCd=0.5;
PCa=2;
scale1=0.15;
scale2=130;


rot=[0*ones(10,1); 10*ones(100,1)];
for i=1:30
    rot=[rot;0;10];
end

N=length(rot);

%base function
sd=15;
x=-pi:pi/100:pi;
u=-pi:pi/100:pi;
PC=circ_vmpdf(x',u,sd);
basevec=[cos(u)' sin(u)'];


%learning
ss_list=[];
w=u*0+w0;
cs_PC=u*0;
for n=1:N-1
    r=rot(n);
    if r>0
        error=0;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    elseif r<0
        error=pi;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    else
        act_CF=0;
    end
    
    cs_PC=cs_PC+PCa.*act_CF;
    
    w=w-r_ltp.*(w-w0)-cs_PC.*r_ltd.*(w-w_lowb).*fCS;
    cs_PC=cs_PC-PCd.*cs_PC;
    ss_PC=w-cs_PC*scale1;
%     w_list(n,:)=w;
%     cs_list(n,:)=cs_PC;
    ss_list(n,:)=ss_PC;
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

% figure; hold on
% subplot(3,2,1);hold on;
% bs_circle=ss_list(end,:)'.*basevec;
% % plot(bs_circle(:,2),bs_circle(:,1),'-k','color',[0.5,0.8,0.7],'linewidth',2);
% fill(sin(x),cos(x),[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.1);
% for i=1:length(u)-1
%     plot(bs_circle(i:i+1,2),bs_circle(i:i+1,1),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',3,'markersize',20);
% end
% axis([-2 2 -2 2],'square')
% 
% subplot(3,2,3);hold on;
% plot(hand*scale2,'-k');
% axis([0 N -5 25])

figure
subplot(2,2,1);hold on;
%label
plot(-cs_list,'-k','color',[0.9,0.75,0],'linewidth',2,'markersize',6);
%stable
plot(w_list,'-k','color',[0.3,0.5,0.9],'linewidth',2,'markersize',6);
%hand
plot(hand*scale2,'-k','color',[0.4,0.6,0.6],'linewidth',2,'markersize',6);

%%
%Exp 3 Antergrade interference
%initial pramater

close all
clear

%initial pramater
error_d=0;
w0=1;
w_lowb=0.15;
r_ltd=0.1;
r_ltp=0.018;
fCS=0.5;
PCd=0.5;
PCa=2;
scale1=0.15;
scale2=130;


rot=[0*ones(10,1); 10*ones(100,1);-10*ones(200,1)];
N=length(rot);

%base function
sd=15;
x=-pi:pi/100:pi;
u=-pi:pi/100:pi;
PC=circ_vmpdf(x',u,sd);
basevec=[cos(u)' sin(u)'];


%learning
ss_list=[];
w=u*0+w0;
cs_PC=u*0;
for n=1:N-1
    r=rot(n);
    if r>0
        error=0;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    elseif r<0
        error=pi;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    else
        act_CF=0;
    end
    
    cs_PC=cs_PC+PCa.*act_CF;
    
    w=w-r_ltp.*(w-w0)-cs_PC.*r_ltd.*(w-w_lowb).*fCS;
    cs_PC=cs_PC-PCd.*cs_PC;
    ss_PC=w-cs_PC*scale1;
    w_list(n,:)=w;
    cs_list(n,:)=cs_PC;
    ss_list(n,:)=ss_PC;
    DCN=mean(ss_PC'.*basevec);
    [dtheta, drho]=cart2pol(DCN(1),DCN(2));
    hand(n)=-drho.*cos(dtheta);
    
    
end

figure; hold on

subplot(3,2,1);hold on;
plot(hand*scale2,'-k');
axis([0 N -30 30])

%%
%Exp 5
%initial pramater

close all
clear

%initial pramater
error_d=0;
w0=1;
w_lowb=0.15;
r_ltd=0.1;
r_ltp=0.018;
fCS=0.5;
PCd=0.5;
PCa=2;
scale1=0.15;
scale2=130;
scale3=0.2;
tim=1:10;
yt=0*[tim tim tim tim tim tim tim tim];

rot=[0*ones(10,1); 10*ones(100,1);-10*ones(12,1);0*ones(30,1)];

r=[10 10 10 10 10 10 10 10 -10 -10 -10 -10 -10 -10 -10 -10 ];
for i=1:20
    rot=[rot; r'];
end

rot=[rot;0*ones(10,1); 10*ones(200,1)];

N=length(rot);

%base function
sd=15;
x=-pi:pi/100:pi;
u=-pi:pi/100:pi;
PC=circ_vmpdf(x',u,sd);
basevec=[cos(u)' sin(u)'];

% plot(x,y)
figure; hold on
set(0,'defaultfigurecolor','w')
%learning
ss_list=[];
w=u*0+w0;
cs_PC=u*0;
for n=1:N-1
    r=rot(n);
    if r>0
        error=0;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    elseif r<0
        error=pi;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    else
        act_CF=0;
    end
    
    cs_PC=cs_PC+PCa.*act_CF;
    
    w=w-r_ltp.*(w-w0)-cs_PC.*r_ltd.*(w-w_lowb).*fCS;
    cs_PC=cs_PC-PCd.*cs_PC;
    ss_PC=w-cs_PC*scale1;
    w_list(n,:)=w;
    cs_list(n,:)=cs_PC;
    ss_list(n,:)=ss_PC;
    DCN=mean(ss_PC'.*basevec);
    [dtheta, drho]=cart2pol(DCN(1),DCN(2));
    hand(n)=-drho.*cos(dtheta);
    
    
end

figure; hold on

%PC in polar
subplot(2,2,1);hold on;
bs_circle=ss_list(end,:)'.*basevec;
% plot(bs_circle(:,2),bs_circle(:,1),'-k','color',[0.5,0.8,0.7],'linewidth',2);
fill(sin(x),cos(x),[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.1);
for i=1:length(u)-1
    plot(bs_circle(i:i+1,2),bs_circle(i:i+1,1),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',3,'markersize',20);
end
axis([-2 2 -2 2],'square')

%PC
subplot(2,2,2);hold on;
bfs=[ss_list(end,:).*circ_vmpdf(x',u,sd)]';
for i=1:10:length(u)
    plot(x,bfs(i,:),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
end

axis([-pi pi 0 2])



% learning curve
subplot(2,2,3);hold on;
plot(hand*scale2,'.k');
axis([0 N -30 30])
% learning vs relearning
subplot(2,2,4);hold on;
plot(hand(end-200:end)*scale2,'-k','color',[0.5,0.8,0.7],'linewidth',2);
plot(hand(10:110)*scale2,'-k','color',[0.4,0.6,0.6],'linewidth',2);
axis([0 100 -5 30])




%%
%TBT Design

close all
clear

%initial pramater
error_d=0;
w0=1;
w_lowb=0.15;
r_ltd=0.1;
r_ltp=0.018;
fCS=0.5;
PCd=0.5;
PCa=2;
scale1=0.15;
scale2=130;
scale3=0.2;



tim=1:10;

yt=0*[tim tim tim tim tim tim tim tim];

rot=[]%[0*ones(10,1); 10*ones(100,1);-10*ones(12,1);0*ones(28,1); 10*ones(100,1)];

r=[10 -10 10 -10 ];
for i=1:300
    rot=[rot; r(randperm(4))'];
end

%rot=[rot;0*ones(10,1); 10*ones(200,1)];

N=length(rot);

%base function
sd=15;
x=-pi:pi/100:pi;
u=-pi:pi/100:pi;
PC=circ_vmpdf(x',u,sd);
basevec=[cos(u)' sin(u)'];


%learning
ss_list=[];
w=u*0+w0;
cs_PC=u*0;
for n=1:N-1
    r=rot(n);
    if r>0
        error=0;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    elseif r<0
        error=pi;
        act_CF=abs(circ_vmpdf(error,u,sd)-error_d)./max(circ_vmpdf(0,u,sd));
    else
        act_CF=0;
    end
    
    cs_PC=cs_PC+PCa.*act_CF;
    w=w-r_ltp.*(w-w0)-cs_PC.*r_ltd.*(w-w_lowb).*fCS;
    cs_PC=cs_PC-PCd.*cs_PC;
    
    ss_PC=w-cs_PC*scale1;
    
    w_list(n,:)=w;
    cs_list(n,:)=cs_PC;
    ss_list(n,:)=ss_PC;
    
    DCN=mean(ss_PC'.*basevec);
    [dtheta, drho]=cart2pol(DCN(1),DCN(2));
    hand(n)=-drho.*cos(dtheta);
    
    
end


% subplot(2,2,1);hold on;
% bs_circle=ss_list(end,:)'.*basevec;
% % plot(bs_circle(:,2),bs_circle(:,1),'-k','color',[0.5,0.8,0.7],'linewidth',2);
% fill(sin(x),cos(x),[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.1);
% for i=1:length(u)-1
%     plot(bs_circle(i:i+1,2),bs_circle(i:i+1,1),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',3,'markersize',20);
% end
% axis([-2 2 -2 2],'square')
%
% subplot(2,2,3);hold on;
% plot(hand*scale2,'.k');
% axis([0 N -30 30])

figure(115)
m=0;
for n=[1 2 3 4 50 51 53 54]
    m=m+1;
    if m<5
        subplot(2,2,1);hold on;
    else
        subplot(2,2,2);hold on;
    end
    bfs=[1-w_list(n,:)]';
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
m=0;
for n=[1 2 3 4 51 52 53 54]
    m=m+1;
    if m<5
        subplot(2,2,1);hold on;
    else
        subplot(2,2,2);hold on;
    end
    bfs=[2-cs_list(n,:)]';
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


% figure; hold on
% subplot(2,2,4);hold on;
% plot(hand(end-100:end)*scale2,'-k','color',[0.5,0.8,0.7],'linewidth',2);
% plot(hand(10:110)*scale2,'-k','color',[0.4,0.6,0.6],'linewidth',2);
% axis([0 100 -5 30])
% 
% % subplot(2,2,1);hold on;
% % bfs=[ss_list(100,:).*circ_vmpdf(x',u,sd)]';
% % for i=1:length(u)
% % plot(x,bfs(i,:),'-r','color',[0.5+0.5*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,0.9],'linewidth',2);
% % end
% %
% % axis([-pi pi 0 2])
% 
% subplot(2,2,2);hold on;
% bfs=[ss_list(end,:).*circ_vmpdf(x',u,sd)]';
% for i=1:10:length(u)
%     plot(x,bfs(i,:),'-r','color',[0.6+0.4*abs(u(i))/pi,0.9-0.3*abs(u(i))/pi,1-0.15*abs(u(i))/pi],'linewidth',1.5);
% end
% 
% axis([-pi pi 0 2])

dccw=find(rot>0);
dcw=find(rot<0);
hand=hand*scale2;
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
axis([0 2 0 5])

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