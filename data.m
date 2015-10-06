function [newrain,newPET,newlen,newQ,daytempmean,newEV,dlhour,awv]=data(newscale,day_start,day_end,tt,metdata)
%metdata is the matrix with all the data over the whole time period
%load metdata.csv; % <-- moved loading to main code
%re-scale from original scale "orscale" to newscale "newscale". 
%scales are in minutes: 1 day=1440 mins

%This is only the time interval I'm intersted in.
date=datenum(metdata(:,1),metdata(:,2),metdata(:,3));
metdata=metdata(find(date>=day_start & date<=day_end),:);

orlen=length(metdata);
orscale=5; %[min]
daymin=1440; %number of minutes in a day
daystep=daymin/orscale;
lenstep=newscale/orscale; % in terms of the original scale
newlen=orlen*orscale/newscale;

% metdata is has a 5mins resolution
year=metdata(:,1);
month=metdata(:,2);
day=metdata(:,3);
rain=metdata(:,4); %[mm]
temp=metdata(:,5); %[C]
press=metdata(:,6); %[mbar]
solrad=metdata(:,7); %[KJ m^(-2)]
hum=metdata(:,8); %[%]
wind=metdata(:,9); %[m/s]
Qdata=metdata(:,10); %[m^3]
%daily variables
dayyear=zeros(1,orlen*orscale/newscale);
daymonth=zeros(1,orlen*orscale/newscale);
dayday=zeros(1,orlen*orscale/newscale);

daytempmean=zeros(1,orlen*orscale/newscale);
daytempmax=zeros(1,orlen*orscale/newscale);
daytempmin=zeros(1,orlen*orscale/newscale)+1000;
daypressmean=zeros(1,orlen*orscale/newscale);
daysolradmean=zeros(1,orlen*orscale/newscale);
dayhummean=zeros(1,orlen*orscale/newscale);
daywindmean=zeros(1,orlen*orscale/newscale);

% EXTRACT DAILY DATA

%daily max and min temp
for i = 1:orlen*orscale/newscale
    for j=i*lenstep-lenstep+1:i*lenstep
        dayyear(i)=dayyear(i)+year(j);
        daymonth(i)=daymonth(i)+month(j);
        dayday(i)=dayday(i)+day(j);
        daytempmean(i)=daytempmean(i)+temp(j);
        daytempmax(i)=max(daytempmax(i),temp(j));
        daytempmin(i)=min(daytempmin(i),temp(j));
        daypressmean(i)=daypressmean(i)+press(j);
        daysolradmean(i)=daysolradmean(i)+solrad(j);
        dayhummean(i)=dayhummean(i)+hum(j);
        daywindmean(i)=daywindmean(i)+wind(j);
    end
        dayyear(i)=dayyear(i)/lenstep;
    daymonth(i)=daymonth(i)/lenstep;
    dayday(i)=dayday(i)/lenstep;
    daytempmean(i)=daytempmean(i)/lenstep;
   % daytempmean(i)=(daytempmax(i)+daytempmin(i))/2;
    daypressmean(i)=daypressmean(i)/lenstep;
    dayhummean(i)=dayhummean(i)/lenstep;
    daywindmean(i)=daywindmean(i)/lenstep;
end

% COMPUTE PET WITH PENMNAN-MONTEITH
% EV is the evaporation from open water computed with Shuttleworth (1993)
% Convert the variables that need to be converted:
daypressmean1=daypressmean*0.1; %[kPa]
daysolradmean1=daysolradmean/1000;
% PET and EV are in mm/day
%dlhour is the number of day light hours and awv is the actual water
%vapor in kPa

%Function that calculates using penmanmonteith.
[PET,EV,dlhour,awv]=PenmanMonteith(daytempmean,daytempmax,daytempmin,daypressmean1,dayhummean,daywindmean,newscale,dayday,daymonth,dayyear,tt);

% Transform to the new time resolution "newscale"
newrain=zeros(1,newlen);
newQ=zeros(1,newlen);
tmp=1:newlen;
tmp1=0;
for i = 1:newlen
    for j=tmp(i)*lenstep-lenstep+1:tmp(i)*lenstep
        newrain(i)=newrain(i)+rain(j);
        newQ(i)=newQ(i)+Qdata(j);
    end
end


%Transform PET and EV from daily to newscale resolution
% y=ones(1,daymin/newscale);
% s=PET'*y;
% s1=EV'*y;
% s=s/(daymin/newscale);
% s1=s1/(daymin/newscale);
% s=s';
% s1=s1';
% newPET=reshape(s,1,(daymin/newscale)*length(PET));
% newEV=reshape(s1,1,(daymin/newscale)*length(EV));
newPET=PET;
newEV=EV;






end


