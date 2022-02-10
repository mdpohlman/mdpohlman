%Mason Pohlman, 11-12-21, HORIZON Group
%FLDI Processor, used for a 2-pt FLDI system (has velocity calculation)
%Help from James Chism

%% Quick Look FLDI
clc; clear variables; close all;

%% File Path
path = 'G:\Mach 4\WavyWall\FLDI\November_2021\Data\';
% % % Photodiode Data % % %
filename = 'Run_5.mat';
File = [path,filename];

%% Tunnel Conditions
% % % Mach Number % % %
M = 4;
% % % Stagnation Pressure % % %
P0Psi = 69*.93;                         %(psia)
% % % Stagnation Temperature % % %
T0_K = 300;                         %(deg K)
% % % Tunnel State Calculator % % %
[PstaticPa,PstaticPsi,Uinf,Re_m] = TunnelState(M,T0_K,P0Psi);
% % % Characteristic Frequency % % %
Del99 = 4.0/1000;                   %(m)
CharFreq = Uinf/Del99;

%% Setup Dimensions
% % % Beam Spacing % % %
BeamSep = 1.09e-3;                  %(m)
% % % Beam Diameter % % %
BeamDiam = 80e-3;                   %(m)
% % % Sample Rate % % %
fs = 2e6;
% % % Number of Samples % % %
Rec_Time = 1;
Samples = floor(Rec_Time*fs);

%% Freestream Data
FS_filename = 'Run_1.mat';
FS_file = [path,FS_filename];
FS_Data = importdata(FS_file);
FS_MS = FS_Data.LogFLDI(160000:400000,:)-mean(FS_Data.LogFLDI(160000:400000,:));

%% Temporal Data
% % % Import Data % % %
Data = importdata(File);
MSData = [Data.LogFLDI-mean(Data.LogFLDI),Data.LogFLDI_1-mean(Data.LogFLDI_1),...
    Data.LogFLDI_2-mean(Data.LogFLDI_2),Data.LogFLDI_3-mean(Data.LogFLDI_3),...
    Data.LogFLDI_4-mean(Data.LogFLDI_4),Data.LogFLDI_5-mean(Data.LogFLDI_5)];
MSData = highpass(MSData(160000:400000,:),1e4,fs);
time = ((1:length(MSData(:,1)))./fs).*1000;   %(ms)

%% Signal Check
figure()
plot(time,MSData(:,1),'linewidth',1.8)
hold on
plot(time,MSData(:,6),'linewidth',1.8)
hold off
title('Raw Signal','Fontsize',16,'fontweight','bold')
xlabel('Time [ms]','Fontsize',14)
ylabel('Voltage [V]','Fontsize',14)
set(gcf,'color','white')
box on

%% Frequency Data
% % % Analysis of Photodiode Data % % %
[PSD,F,Gxxfs,Ffs,Gxx,dF,dt,Velocity,Lags,Coef] = DiodeAnalysis(MSData,FS_MS,Samples,fs,BeamSep);
% % % Signal Velocity Normalized by Ideal FS Velocity % % %
V_Tunn = Velocity/641;
% % % Signal Velocity Normalized by Measured FS Velocity % % %
% V_FLDI = Velocity/FLDI_Uinf;

%% Plots
FigureGeneration(F,Gxx,Gxxfs,Ffs,CharFreq);

%% Functions
% % % Tunnel State % % %
function [PstaticPa,PstaticPsi,Uinf,Re_m] = TunnelState(M,T0_K,P0Psi)
% James Chism
% Reynolds number is calculated per meter. To convert it, the user must
% input their characteristic length in meters and inches.

%% Change Log
% % % 8/20/2021 % % %
% Matched with Mark's calculator, Mu equation was wrong.

%% Input Values
% % Tunnel mach number
MP2 = M^2;

% % Stagnation pressure. For Mach 2 use 35 psi, for Mach 4 use 43.4
% P_Input = input('Stagnation Pressure (psi) = ');
% P0Psi = 41;
P0Pa = P0Psi*6894.76;

% % Specific Heat 
Gamma = 1.4;
GM1 = Gamma-1;
GP1 = Gamma+1;

% % Gas Constant (J/kg*K)
R = 287.058;

%% Calculated Values
% % % Total Pressure % % %
P02_P01 = ((GP1*MP2)/(GM1*MP2+2))^(Gamma/GM1)*(GP1/(2*Gamma*MP2-GM1))^(1/GM1);

% % % Area Ratio % % %
A_Astar = ((GP1/2)^(-GP1/(2*GM1)))*(((1+((GM1/2)*MP2))^(GP1/(2*GM1)))/M);

% % % Total Density % % %
Rho0 = P0Pa/(R*T0_K);     %(kgm^3)

% % % Static Pressure % % %
PstaticPa = ((1+(GM1/2)*MP2)^(-Gamma/GM1))*P0Pa;
PstaticPsi = ((1+(GM1/2)*MP2)^(-Gamma/GM1))*P0Psi;

% % % Static Density % % %
Rho = ((1+(GM1/2)*MP2)^(-1/GM1))*Rho0;       %(kg/m^3)

% % % Static Temperature % % %
T = ((1+(GM1/2)*MP2)^-1)*T0_K;

% % % Dynamic Viscosity % % %
Mu = (0.00001716*((T/273.15)^1.5)*(273.15+110.4))/(T+110.4);         %(Pa/s)

% % % Kinematic Viscosity % % %
Nu = Mu/Rho;        

% % % Speed of Sound % % %
a = sqrt(Gamma*PstaticPa/Rho);

% % % Velocity % % %
Uinf = a*M;

% % % Unit Reynolds Number % % %
Re_m = (Rho*Uinf)/Mu;       %(m^-1)
% Re_ft = Re_m*(Char_Length_m/Char_Length_in)*12;
end
% % % Photodiode Analysis % % %
function [PSD,F,Gxxfs,Ffs,Gxx,dF,dt,Velocity,Lags,Coef] = DiodeAnalysis(MSData,FS_MS,Samples,fs,BeamSep)
Win= hanning(floor(Samples/100));

% % % Freestream % % %
[PSDfs,Ffs] = pwelch(FS_MS,Win,[],[],fs);
dFfs = Ffs(2)-Ffs(1);
Gxxfs = PSDfs.*dFfs./var(FS_MS);

% % % Model Run % % %
[PSD,F] = pwelch(MSData,Win,[],[],fs);

dF = F(2)-F(1);

Gxx = PSD.*dF./var(MSData);

%% Diode Velocity
% % % Model % % %
[Coef,Lags] = xcorr(MSData(:,1),MSData(:,6),'coef');
[~,loc] = max(abs(Coef));
dt = abs(Lags(loc)/fs);
Velocity = (5*BeamSep)/dt;

end
% % % Figure Generation % % %
function FigureGeneration(F,Gxx,Gxxfs,Ffs,CharFreq)
% % % Frequency Subplot % % %
fig = figure();
% % % Ch 1 % % %
subplot(3,2,1)
loglog(F(18:end),Gxx(18:end,1),'k','Linewidth',1.0)
axis([F(18),1e6,1e-6,1e-2])
text(350,3e-6,'Ch 1','HorizontalAlignment','right')
box on

% % % Ch 2 % % %
subplot(3,2,2)
loglog(F(18:end),Gxx(18:end,2),'k','Linewidth',1.0)
axis([F(18),1e6,1e-6,1e-2])
text(350,3e-6,'Ch 2','HorizontalAlignment','right')
box on

% % % Ch 3 % % %
subplot(3,2,3)
loglog(F(18:end),Gxx(18:end,3),'k','Linewidth',1.0)
axis([F(18),1e6,1e-6,1e-2])
text(350,3e-6,'Ch 3','HorizontalAlignment','right')
box on

% % % Ch 4 % % %
subplot(3,2,4)
loglog(F(18:end),Gxx(18:end,4),'k','Linewidth',1.0)
axis([F(18),1e6,1e-6,1e-2])
text(350,3e-6,'Ch 4','HorizontalAlignment','right')
box on

% % % Ch 5 % % %
subplot(3,2,5)
loglog(F(18:end),Gxx(18:end,5),'k','Linewidth',1.0)
axis([F(18),1e6,1e-6,1e-2])
text(350,3e-6,'Ch 5','HorizontalAlignment','right')
box on
% % % Ch 6 % % %
subplot(3,2,6)
loglog(F(18:end),Gxx(18:end,6),'k','Linewidth',1.0)
axis([F(18),1e6,1e-6,1e-2])
text(350,3e-6,'Ch 6','HorizontalAlignment','right')
box on

% % % Plot Options % % %
this=axes(fig,'visible','off'); 
this.Title.Visible='on';
this.XLabel.Visible='on';
this.YLabel.Visible='on';
ylabel(this,'G(f)\cdotdf/\sigma^2','Fontsize',14);
xlabel(this,'Frequency [Hz]','Fontsize',14);
title(this,'Freestream Spectra','Fontsize',16,'fontweight','bold');
set(gcf,'color','white')

% % % Frequency Comparison % % %
figure()
% % % Ch 1 % % %
loglog(F(18:end),Gxx(18:end,1)*10e3,'k','Linewidth',1.0)
hold on
% % % Ch 2 % % %
loglog(F(18:end),Gxx(18:end,2)*10e2,'b','Linewidth',1.0)
% % % Ch 3 % % %
loglog(F(18:end),Gxx(18:end,3)*10e1,'g','Linewidth',1.0)
% % % Ch 4 % % %
loglog(F(18:end),Gxx(18:end,4)*10e0,'y','Linewidth',1.0)
% % % Ch 5 % % %
loglog(F(18:end),Gxx(18:end,5)*10e-1,'r','Linewidth',1.0)
% % % Ch 6 % % %
loglog(F(18:end),Gxx(18:end,6)*10e-2,'m','Linewidth',1.0)
% % % Characteristic Frequency % % %
LineText = ['F_c = ',num2str(round(CharFreq/1000,0)),' [kHz]'];
xline(CharFreq,'--k',LineText,'LabelHorizontalAlignment','left',...
    'LabelVerticalAlignment','bottom','LabelOrientation','Horizontal','fontsize',11);
hold off
legend('Ch 1','Ch 2','Ch 3','Ch 4','Ch 5','Ch 6','location','southwest','NumColumns',2,'fontsize',11)
ylabel('G(f)\cdotdf/\sigma^2','Fontsize',14);
xlabel('Frequency [Hz]','Fontsize',14);
title('Single Plot Comparison','Fontsize',16,'fontweight','bold')
axis([F(18),1e6,10e-11,10e0])
box on
grid on
set(gcf,'color','white')

% % % Freestream Comparison % % %
figure()
% % % Freestream % % %
loglog(Ffs(18:end),Gxxfs(18:end,1),'k','Linewidth',1.0)
hold on
% % % Model % % %
loglog(F(18:end),Gxx(18:end,1),'r','Linewidth',1.8)
% % % Characteristic Frequency % % %
LineText = ['F_c = ',num2str(round(CharFreq/1000,0)),' [kHz]'];
xline(CharFreq,'--k',LineText,'LabelHorizontalAlignment','left',...
    'LabelVerticalAlignment','bottom','LabelOrientation','Horizontal','fontsize',11);
hold off
legend('Freestream','Model','location','southwest','fontsize',11)
ylabel('G(f)\cdotdf/\sigma^2','Fontsize',14);
xlabel('Frequency [Hz]','Fontsize',14);
title('Freestream Model Comparison','Fontsize',16,'fontweight','bold')
axis([F(18),1e6,10e-11,10e0])
set(gcf,'color','white')
box on
grid on

end
