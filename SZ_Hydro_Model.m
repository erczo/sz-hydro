% SZ HYDRO REW MODEL 
% Author: Stefano Zenardo
% Date: 2013-04-05
% 
% Modified by:  Collin Bode and Daniella Rempe
% Date: 2013-
%
% Description:  
% Simple hydrology model that uses Relative Elementary Watersheds (REW) as
% its fundamental hydrologic unit.  Each REW consists of a polygon and a
% line representing the unit and the stream channel it drains into.
% 
% WATER BALANCE
% 1-dimensional 3-layer model describes water balance within the REW.
% Layer 1. "root zone" = initial storage + rainfall - percolation to Layer
% 2.
% Layer 2. "fractured rock zone" = initial storage + percolation from Layer
% 1 - outflow to stream - percolation to Layer 3.  This is rapid pulsed
% flow.
% Layer 3. "bedrock zone" = initial storage + percolation from Layer 2 -
% outflow to stream.  This is the slow baseflow.
% 
% KINEMATIC WAVE 
% Water movement in streams is done using kinematic wave.  Stream is
% segmented into reaches.  Watershed is processed from upper elevations on down based on FROM TO
% discharge dependencies.  
%
% Note1: Stream _must_ be ordered by REW's at the top of
% the watershed flowing downstream.  Will fail if ordering is incorrect.
%
% Note2: Stream movement timestep must be short enough that the water
% cannot move through the reach within one timestep.  Equations will go
% unstable if that is the case.
%
% DEPEDENCIES
%
% Matlab Codes:
% SZ_Hydro_Model.m main code, this file
%
% Functions:
% net_analysis.m   performs stream channel movement
% data.m           interpolates input meteorlogical values to timestep 
% PenmanMonteith.m
% rmse.m
% r2.m
% 
% Data Files:
% cal_TMP.mat
% datamat.csv
% elder30mtree150.dat
% net150.csv
%---------------------------

close all
time=cputime;

% Load all data files
 
% Meteorological Data
% datamat is the matrix with all the meteorological and discharge data over the whole time period
load datamat.csv;
metdata=datamat;

% Stream Network data
% net1 =  
% tree = 
%Make sure the input is ordered properly
%for the Elder use flag1=1, for the Eel river use flag1=0
flag1=1;
if flag1
    csvread net150.csv;
    net1 = ans;
    tree = load('elder30mtree150.dat');
else
    csvread Eelnet.csv;
    net1 = ans;
    tree = load('Eeltree.dat');
end


%-----------------------------

if exist('rain')	% Do not load the cal_TMP or Meteorological datasets if they are already loaded.
else
    %load calibrationNetwork;
    %x=cal(find(cal(:,9)>0.9 & cal(:,10)>0.96 & cal(:,11)<1.67),:);

	instchan=[];
	area=[];
	% Calibration matrix
    load cal_TMP;
    cal=cal_TMP;
    % cal -> [ 1n 2Z 3s1 4sw 5s0 6log10(Ksatday) 7kday 8kgwday 9c 10f 11b 12cw
    %          13m 14Ens 15R2 16Rmse 17CRM 18max(d_daily)]

	day_start=datenum(2010,9,1);
	day_end=datenum(2011,9,1);
	Nday=day_end-day_start+1;   %number of day of simulation
	A=16.87; %[km^2]
	dt=15; %[min] Time step for the hillslope
	dtk=1; %[min] new time step for the kinematic wave model in the river network
	dt_day=1440/dtk;
	dt_daydata=1440/dt;
	% Questo e' complicato, tt e' il vettore dei punti medi negli intervalli di
	% tempo considerati serve per calcolare PET.
	%This is the vector of the middle points in the time intervals
	%considered to calculate the potential evapotranspiration.
	tt=zeros(1,24*60/dt);
	tt(1)=-dt/60/2;
	for i=2:length(tt)
		tt(i)=tt(i-1)+dt/60;
	end
	tt=tt+60/dt;
	tmp=ones(Nday+1,1)*tt;
	tt=reshape(tmp',1,length(tt)*(Nday+1));

	%This function takes the data and interoplates or averages them in
	%order for them to have the same time-step.
	[rain,PET,N,Qdata,Ta,EV,dlhour,awv]=data(dt,day_start,day_end,tt,metdata);
	%dlhour is the number of day light hours and avp is the actual water
	%vapor in kPa
	%N=number of time steps
	%Qdata is in m^3 need to convert it in mm:
	Qdata=Qdata/A/1000;
	T=dt*N; %time of simulation [min]
	%PET and EV are in mm/day, I need to transform it into mm in the
	%interval dt:

	PET=PET*dt/1440;
	EV=EV*dt/1440;
	%Transform water vapor from kPa to Pa:
	awv=awv*1000;
	%Transform dlhour at the daily scale: daydlhour
	tmp=1:length(dlhour);
	daydlhour=interp1(tmp,dlhour,1:length(dlhour)/(day_end-day_start):length(dlhour));
end

x=find(cal(:,14)>0.80 & cal(:,15)>0.9 & cal(:,17)<0.05 & cal(:,17)>-0.05);
%x=find(cal_Sep2012(:,15)>0.9);

y=1;
Sx=0;
Qx=0;
Vx=0;
Jex=0;
Qinx=0;
Qtotx=0;
ETx=0;
rainx=0;
snZx=0;

%CALL THE FUNCTION FOR THE STREAM NETWORK GEOMETRY AND CONNECTIVITY
%Channel width is computed using Montgomery and Gran 2001 WRR (Table 1),
%w=c*A^b, w is in meters and A is in square meters.
c=cal(x(y),12);%0.01;
b=cal(x(y),11);%0.39;
[link_no,downstream_link_no,upstream_link_no1,upstream_link_no2,order,len,...
    magnitude,downstream_contr_area,slope,upstream_contr_area...
    ,watershed_no,watershed_area,wid,new_order]= net_analysis(c,b,net1,tree);

%ii=1;
%----- parameters for soil moisture dynamics ------------------------------
n=cal(x(y),1);%x(ii,1); %0.1;      %porosity
Z=cal(x(y),2);%3000;    %[mm]
s1=cal(x(y),3);%x(ii,3);    %above s1 plant evaporates at potential ET
sw=cal(x(y),4);%x(ii,4);    %wilting point
s0=cal(x(y),5);%x(ii,5);    % initial soil water content
Ksat=10^cal(x(y),6);%x(ii,6); %0.5 %[mm/min] hydraulic conductivity at saturation
%k =1/(30*24*60);%/sqrt(area/sum(watershed_area)); %1/(x(ii,7)*24*60);  %subflow rate [1/min] (parameter of the exponential IUH, mean IUH equal to 1/k)
k =1/(cal(x(y),7)*24*60);%/sqrt(area/sum(watershed_area)); %1/(x(ii,7)*24*60);  %subflow rate [1/min] (parameter of the exponential IUH, mean IUH equal to 1/k)
kgw=1/(cal(x(y),8)*24*60);%/sqrt(area/sum(watershed_area));   %1/(x(ii,8)*24*60); %[1/min]
c=cal(x(y),9);%x(ii,9);      %3.4801;     %Clapp exponent
Recharge=6.7294e-004; %[mm/min]
f=cal(x(y),10);%0.06;       %x(ii,10);
% Change of dimensions

Ksat=Ksat*dt;
k=k*dt;
kgw=kgw*dt;
Recharge=Recharge*dt;






%**************************************************************************
% Numerical Simulation --> Soil moisture dynamics
%**************************************************************************
%NB! no surface runoff is taken in account

%----- fluxes
Je=zeros(1,N);  %percolation (as in the paper)
ET=zeros(1,N);  %evapotranspiration
Qsub=zeros(1,N);   %subsurface runoff
R=zeros(1,N); %groundwater recharge
Qgw=zeros(1,N); %groundwater
Qtotal=zeros(1,N); %total discharge per unit area



%----- storage
s=zeros(1,N);   %soil moisture [-]
ST=zeros(1,N);  %storage nel subsurface [m3/m2]
Sgw=zeros(1,N); %groundwater storage [m3/m2]

%Initial conditions
s(1)=s0;       %initial soil moisture value


%----- Numerical Simulation -----------------------------------------------
%   rain=rain/dt; %[mm/min]
%   PET=PET/dt  ; %[mm/min]


for t=1:N  %for loop on the number of timesteps
    
    %upper storage
    ET(t)=max(0,min(PET(t),PET(t)*((s(t)-sw)/(s1-sw)))); %[dt*mm/min]
    %r=[s1 s0 sw Ksat k kgw mean(PET) n Z]
    % keyboard
    Je(t)=Ksat*s(t)^c; %[dt*mm/min]
    
    s(t+1)=min(1,s(t)+(rain(t)-ET(t)-Je(t))/(n*Z));  %soil moisture cannot exceed one
    %R(t)=f*Je(t);  %percolation to the third layer
    R(t)=min(ST(t),Recharge); %fixed deep percolation rate from second to third layer
    
    Qsub(t)=k*ST(t); %[dt*mm/min]=[dt*m3/m2/min]
    
    % Water balance in the second layer
    ST(t+1)=ST(t)+(Je(t)-Qsub(t)-R(t));
    
    %Third layer
    Qgw(t)=kgw*Sgw(t);
    Sgw(t+1)=Sgw(t)+(R(t)-Qgw(t));
    
    
    dl=1; %uni
    
    Qtotal(t)=(Qsub(t)+Qgw(t))*0.001*dl; %[m] 0.001 is to go from mm to m.
    %The variables with x are only used to do the mass balance check.
    Jex=Jex+Je(t)*0.001;   %[m]
    %Qinx=Qinx+(Qsub(t)+Qgw(t))*0.001; %[m]
    Qtotx=Qtotx+Qtotal(t); %[m]
    ETx=ETx+ET(t)*0.001; %[m]
    rainx=rainx+rain(t)*0.001; %[m]
end

snZx=snZx+s(N)*n*Z*0.001; %[m^3]
Sx=Sx+(ST(N)+Sgw(N))*0.001; %[m^3]



%**************************************************************************
% KINEMATIC WAVE MODEL
% REDUCE THE TIME STEP TO AVOID INSTABILITY IN THE KINEMATIC WAVE MODEL
% See the time-step dtk at the top of the page.
% dt/dtk MUST BE INTEGER
%**************************************************************************

TMP=zeros(dt/dtk,N);
for i=1:dt/dtk;
    TMP(i,:)=Qtotal(:)/(dt/dtk);
end
Qtotal1=reshape(TMP,N*(dt/dtk),1);
dt_dayk=1440/dtk;
Nk=length(Qtotal1);

Qtot=zeros(1,Nk); %total discharge
d=zeros(length(watershed_area),Nk+1);
V=zeros(length(watershed_area),Nk+1);
Qin_net=zeros(length(watershed_area),Nk+1);
Qout_net=zeros(length(watershed_area),Nk+1);
Hill_discharge=zeros(length(watershed_area),Nk+1);

d=d';
V=V';
Qin_net=Qin_net';
Qout_net=Qout_net';
Hill_discharge=Hill_discharge';

theta=0; %For theta=0 the cross-section is rectangular
% if theta == 0;
%     logictheta=1;
% else
%     logictheta=0;
%     'REMEMBER TO CHANGE THE USE OF THETA INSIDE THE CYCLE'
%     pause
% end
TMP=0;
% in is the the index for each link/stream/channel
for in=1:length(watershed_area)
    
    %Parameters Kinematic wave
    
    dl=1; %Not important parameter, keep it equal to one.
    
    m= cal(x(y),13);%0.03; %Manning coefficient s/m^(1/3)
    
   
    m=m/(60*dtk)/(dl^(1/3)); %first transformed in min/m^(1/3) and then transfromed in 1/m^(1/3) by multiplying by the current time scale
    slo=slope(in);
    w=wid(in)*dl;
    L=len(in)*dl;
    area=watershed_area(in)*dl^2;
    downarea=downstream_contr_area(in)*dl^2;
    linkno=link_no(in);
    downlinkno=downstream_link_no(in);
    % In some cases when we have channels of 1 pixels, area or slope might be
    % zero
    if area==0
        area=30*30;
    end
    if slo==0
        slo=0.001
    end
    
    %Initial conditions kinematic wave
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q0=downarea*(f/(kgw*(dtk/dt)))*(Ksat*(dtk/dt)*0.001*dl)*s0^c+downarea*((1-f)/(k*(dtk/dt)))*(Ksat*(dtk/dt)*0.001*dl)*s0^c; %[m^3]
    syms d1;
    F=solve(Q0 - (d1*sqrt(slo)*(w + d1*tan(theta))*((d1*(w +...
        d1*tan(theta)))/(w + 2*d1*sec(theta)))^(2/3))/m,d1);
    %Select only the solution that is real positive. I am sure this can be done
    %in a easier way...
    for iF=1:length(F)
        if isreal(F(iF))
            if double(F(iF))>0
                iFF=iF;
            end
        end
    end
    d(in,1)=double(F(iFF));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Stream channel volume
    V(in,1)=d(in,1)*L*w;
    
    
    Qtot=Qtotal1*area; %[m^3]
    
    
    for t=1:Nk
        Qin_t=Qin_net(t,in)+Qtot(t); %[m^3]
        d_t=d(t,in);
        
        %Qin_net(t,in)=Qin_net(t,in)+Qtot(t); %[m^3]
        
        Hill_discharge(t,in)=Qtot(t); %Probably not necessary
        
        
        Qout_t=(d_t*sqrt(slo)*w*((d_t*w)/(2*d_t + w))^0.67)/(m);
        V_t1=max(0,d_t*L*w+(Qin_t-Qout_t));
        d_t1=V_t1/(L*w); %This is d at t+1
    
        
        % CHECK FOR NUMERICAL INSTABILITY
        lenv=Qout_t/(d_t*w); % Distance run by water in the time interval dt
        %lenv/dt= water velocity
        % Since I'm using an explicit method, if lenv>L we have numerical instability
        tmp=0;
        if lenv>L
            %[lenv/dtk lenv L in t] %lenv/dtk=[m/min]
            %'WARNING!!! NUMERICAL INSTABILITY!!!'
            %lenv/dtk/60 % [m/s]
            %tmp=in
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS PART IS NOT CORRECT THEORETICALLY BUT FIXES THE INSTABILITY
            % PROBLEM AND IT'S OK IF APPLIED TO SHORT CHANNELS. IF INSTABILITY
            % OCCURS IN MANY CHANNELS THE TIME-STEP SHOULD BE REDUCED
            Qout_t=Qin_t;
            V_t1=max(0,d_t*L*w);%+(Qin_net(t,in)-Qout_net(t,in)));
            d_t1=V_t1/(L*w);
            %instchan=[instchan in];
            %if TMP<in
            %    instchan=[instchan in];
            %    TMP=in;
            %end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %Qin_net(t+1,mask)=Qin_net(t+1,mask)+Qout_t;
        
        
        Qout_net(t,in)=Qout_t;
        Qin_net(t,in)=Qin_t;
        d(t,in)=d_t;
        V(t+1,in)=V_t1;
        d(t+1,in)=d_t1;
        
        
        %Qx=Qx+Qout_net(length(watershed_area),t); %[m^3]
    end
    
    % Update the following channels
    if in<length(watershed_area)
        mask=find(link_no==downlinkno);
        tmp3=[0 ;Qout_net(:,in)];
        tmp3(end)=[];
        Qin_net(:,mask)=Qin_net(:,mask)+tmp3;
    end
    
    
    
    
    Qx=sum(Qout_net(:,length(watershed_area))); %[m^3]
    Vx=Vx+V(Nk,in); %[m^3]
    
    in
end
Qinx=Qtotx*sum(watershed_area);
% MASS BALANCE CHECK
if (rainx-snZx-ETx-Jex)/(rainx)>0.001
    'MASS BALANCE NOT SATISFIED IN THE ROOT ZONE'
    pause
end
if (Jex-Sx-Qtotx)/(Jex)>0.001
    'MASS BALANCE NOT SATISFIED IN THE TRANSPORT ZONE'
    pause
end
if (Qinx-Vx-Qx)/(Qinx)>0.001
    'MASS BALANCE NOT SATISFIED IN THE RIVER NETWORK ZONE'
    pause
end

elapsed=cputime-time

%Convert back to the daily time-scale
Qout_net(1,:)=[];
Qin_net(1,:)=[];
Hill_discharge(1,:)=[];
V(1,:)=[];
d(1,:)=[];
Qdata_daily=sum(reshape(Qdata,dt_daydata,N/dt_daydata)); %[mm]
Q_daily=sum(reshape(Qout_net(:,length(watershed_area)),dt_dayk,Nk/dt_dayk))*1000/sum(watershed_area); %[mm]
d_daily=mean(reshape(d(:,length(watershed_area)),dt_dayk,Nk/dt_dayk)); %[mm]

%EFFICIENCY INDICATORS
Ens=NashSutcliffe(Qdata_daily,Q_daily)
Rmse=rmse(Qdata_daily,Q_daily)
R2=r2(Qdata_daily,Q_daily)
CRM=(sum(Qdata_daily)-sum(Q_daily))/sum(Q_daily)

if 1
    %PLOTS
    figure
    xData = linspace(day_start,day_end,length(Q_daily));
    subplot(1,3,1)
    plot(Q_daily)
    hold on
    plot(Qdata_daily,'r')
    set(gca,'XTick',xData)
    datetick('x','mmm-yy','keeplimits')
    title(['Ens=',num2str(Ens),'---R2=',num2str(R2)])
    legend('OBSERVED RUNOFF','MODELED RUNOFF','Location','northeast')
    xlabel('DATE');
    ylabel('RUNOFF [mm]');
    subplot(1,3,2)
    plot(d(:,length(watershed_area)))
    set(gca,'XTick',xData)
    datetick('x','mmm-yy','keeplimits')
    title(['y=',num2str(y)])
    legend('MODELED STREAM DEPTH','Location','northeast')
    xlabel('DATE');
    ylabel('DEPTH [m]');
    subplot(1,3,3)
    plot(cumsum(Q_daily))
    hold on
    plot(cumsum(Qdata_daily),'r')
    set(gca,'XTick',xData)
    datetick('x','mmm-yy','keeplimits')
    legend('OBSERVED CUMUL. RUNOFF','MODELED RUNOFF','Location','northeast')
    xlabel('DATE');
    ylabel('CUMULATED RUNOFF [mm]');
    
    %Convert Ta, EV, and awv in the dtk time scale
end

    
if exist('Qout_net.mat','file')
    'SZ_Hydro_Model.m Output files already exist. Exiting.' 
else
    'SZ_Hydro_Model.m Output files do not exist. Exporting.'
    %Convert Ta, EV, and awv in the dtk time scale
    Ta1=repmat(Ta',1,dt/dtk)';
    Ta=Ta1(:);
    awv1=repmat(awv',1,dt/dtk)';
    awv=awv1(:);
    EV1=repmat(EV',1,dt/dtk)';
    EV=EV1(:)/(dt/dtk);
    save('Hill_discharge.mat','Hill_discharge','-MAT');
    save('Qout_net.mat','Qout_net','-mat');
    save('Qin_net.mat','Qin_net','-mat');
    save('V.mat','V','-mat');
    save('d.mat','d','-mat');
    save('Ta.mat','Ta','-mat'); %Atmospheric temperature in C
    save('EV.mat','EV','-mat'); %Evaporation from open water, mm in the time step dt
    save('daydlhour.mat','dlhour','-mat'); %Number of daylight hours at the daily scale
    save('awv.mat','awv','-mat'); %Actual water vapor in Pa
    save('watershed_area.mat','watershed_area','-mat'); % [m^2] sub-catchments area
end
'DONE!'
