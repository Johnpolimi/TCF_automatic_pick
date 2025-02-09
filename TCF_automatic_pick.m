%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TCF automatic pick
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

% Import data: Test_data1 -------------------------------------------------------------
[fileName, pathName] = uigetfile('*.txt','Choose file to import');
if fileName == 0
    return
end
fileNameComplete = [pathName fileName];
[data1]= load(fileNameComplete);
fs=1000;
data2=data1.';   % transpose of matrix

% Set axes ----------------------------------------------------------------
[ny,nx] = size(data2);
dt = 1/fs*1000;
xAx = 1:nx;
timeAx = 0:dt:(ny-1)*dt;

for i=1:nx
    eventw(:,i) = fft(data2(:,i))*dt;
    eventw(1,i) = 0;                                % remove the DC term
    eventf(:,i) = abs(eventw(:,i));
    filter_event(:,i) = ifft(eventw(:,i))/dt;   
end

f = (0:ny-1)*fs/ny;
timexlim=[0 max(timeAx)];
freqxlim=[0 500];

%% Import manual pick data:  Test_data2 -------------------------------------------------------------
[fileName, pathName] = uigetfile('*.txt','Choose file to import');
if fileName == 0
    return
end
fileNameComplete = [pathName fileName];
[pick1]= load(fileNameComplete);

% TCF Development ------------------------------------------------------
sw=0.3;    % the length of the first window [s]
lw=1.2;      % the length of the second and third windows [s]
Nsw = sw*fs;      % Number of samples in the windows
Nlw = lw*fs;

% calculate the E1, E2 and E3; ER12, ER13, CF
E1=zeros(ny,nx);
E2=zeros(ny,nx);
E3=zeros(ny,nx);
ER12=zeros(ny,nx);
ER13=zeros(ny,nx);
CF=zeros(ny,nx);

aj=filter_event.^2;

for ch=1:nx
    for i=Nsw+Nlw:ny-Nsw
        % Calculate average energy
        E1(i,ch)=sum(aj(i:i+Nsw,ch))/Nsw;
        E2(i,ch)=sum(aj(i-Nlw+1:i,ch))/Nlw;
        E3(i,ch)=sum(aj(i-Nlw-Nsw+1:i-Nsw,ch))/Nlw;
        ER12(i,ch)=E1(i,ch)/E2(i,ch);
        ER13(i,ch)=E1(i,ch)/E3(i,ch);
        CF(i,ch)=ER13(i,ch)-ER12(i,ch);   % Not use alpha
    end
end

% linear correction between max(CF) and max(CF)-2sw
F2=CF;
F3=zeros(ny,nx);
pick=zeros(15,2);
for j=1:nx
    [max_DER,I_DER]=max(F2(:,j));
    a1=(F2(I_DER,j)-F2(I_DER-2*Nsw,j))/(2*Nsw);
    b1=F2(I_DER-2*Nsw,j);
    for i=I_DER-2*Nsw:I_DER
        F3(i,j)=F2(i,j)-(a1*(i-(I_DER-2*Nsw))+b1);
    end
    [pick(j,1),pick(j,2)]=min(F3(:,j));
end

pick_time=pick(:,2)-1; % point in pick(2,:), coordinate-1,i.e. 'pick(2,:)-1', is the pick time

%% Figure1
% Using Test_data1: Fracture 1, Channel 1.
for i=1:1
    figure
    subplot(3,1,1)                                 
    plot(timeAx, filter_event(:,i),'k','LineWidth',0.6)
%     grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'off';
%     ax.TickDir = 'out';
    ylabel('Amplitude');
%     legend('X','Location','northeast','Orientation','horizontal')
    xlabel('Time [ms]');
%     xlim(timexlim)
%     xlim(max(abs(xlim)).*[-0.05 1.05])
    xlim([1000 2600])
    ylim(max(abs(ylim)).*[-0.85,0.85])
    hold on
    yline(0,'--','LineWidth',0.6,'color',[0.5 0.5 0.5]);   % 0 amplitude line
    % plot pick4 first, then pick1
%     xline(pick_time(i),'r','LineWidth',1.1)  % plot automatic pick by CF
%     xline(pick1(i),'--b','LineWidth',1.1)  % plot manual pick

    subplot(3,1,2)
    plot(timeAx, CF(:,i),'k','LineWidth',1.1)
    yline(0,'--','LineWidth',0.6,'color',[0.5 0.5 0.5]);   % 0 amplitude line 
%     grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'off';
%     ax.TickDir = 'out';
    ylabel('CF');
%     legend('ER12','ER13','IDER','Location','northeast','Orientation','horizontal')
%     xlim(freqxlim)
    xlabel('Time [ms]');
%     xlim(max(abs(xlim)).*[-0.05 1.05])
    xlim([1000 2600])
    ylim(max(abs(ylim)).*[-0.1,1.1])

    subplot(3,1,3)                                 
    plot(timeAx, F3(:,1),'k','LineWidth',0.7)
    % title(['Fracture1-CH',num2str(1)])
    ylabel('TCF');
    xlabel('Time [ms]');
    % xlim(timexlim)
%     xlim([abs(max(xlim))*(-0.05),abs(max(xlim))*1.05])
    xlim([1000 2600])
    ylim([abs(min(ylim))*(-1.2),max(abs(ylim))*1.2])
    hold on
    yline(0,'--','color',[0.5 0.5 0.5],'LineWidth',0.7);   % 0 amplitude line
    xline(pick_time(i),'r','LineWidth',1.1)
    xline(pick1(i),'--b','LineWidth',1.1)  % plot manual pick
    legend('','','TCF picking','Manual picking','Location','northwest','Orientation','horizontal')
end










