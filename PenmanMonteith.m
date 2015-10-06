function [PETPM,EV,N,ea]=PenmanMonteith(T,Tmax,Tmin,P,RH,u2,newscale,day,month,year,tt)
% COMPUTE PET WITH PENMNAN-MONTEITH USING FAO APPROACH
%T=temperature in [C]
%P=atmospheric pressure [kPa]
%RH=Relative humidity [%]
%u2=wind speed [m s^-1]

%if latitude is in degree: 39° 43', latitude in decimal degree is 39+43/60
latdeg=39.7167; %decimal degrees
latrad=latdeg*pi/180; %radiant
%longitude time zone longtz=120 for pacific zone
longtz=120;
%longitude of the measuremnt site longm
longm=123;
%z=elevation above sea level [m]
z=1000; %[m] PUO ESSERE PIU PRECISO!!!!
% t1 is the length of the calculation period in hour: i.e., 1 for hourly 
% period or 0.5 for a 30-minute period.
t1=newscale/60;
rhow=1000; %[kg/m3]



for i=1:length(T)
    %hour(i)=
    %es=saturation vapor pressure [kPa]
    es(i)=0.61121*exp(T(i)*(18.678-(T(i)/234.5))/(257.14+T(i)));
    %ea=actual vapor pressure [kPa]
    ea(i)=es(i)*RH(i)/100;
    %G=soil heat flux [MJ m^(-2)]. NEGLECTABLE FOR DAILY SCALES
    G=0;
    %Delta=slope vapor pressure curve [kPa C^(-1)]
    Delta(i)=4098*(0.6108*exp(T(i)*17.27/(T(i)+237.3)))/((T(i)+237.3^2));
    %lam=latent heat of vaporization [MJ Kg^(-1)]
    lam=2.45;
    %gam=psychrometric constant [kPa C^(-1)]
    gam(i)=P(i)*1.013*0.001/(lam*0.622);
    % PETPM=Potential evapotranspiration [mm] with Penman-Monteith
    %J=number of the day of the year: 1=January 1
    J(i)=datenum(year(i),month(i),day(i))-datenum(year(i),1,1)+1;
    Gsc=0.0820; %[MJ m^-2 min^-1]
    dr(i)=1+0.033*cos(J(i)*2*pi/365);
    de(i)=0.409*sin((J(i)*2*pi/365)-1.39);
    omegas(i)=acos(-tan(latrad)*tan(de(i)));
    b(i)=2*3.14*(J(i)-81)/364;

    Sc(i)=0.1645*sin(2*b(i))-0.1255*cos(b(i))-0.025*sin(b(i));
    omega(i)=(3.14/12)*((tt(i)+0.06667*(longtz-longm)+Sc(i))-12);
    omega1(i)=omega(i)-(3.14*t1)/24;
    omega2(i)=omega(i)+(3.14*t1)/24;
    %Ra=extraterrestrial radiation [MJ m^-2 day^-1] for daily or longer
    %scale
    Ra(i)=(24*60/pi)*Gsc*dr(i)*((omegas(i))*sin(latrad)*sin(de(i))+cos(latrad)*cos(de(i))*sin(omegas(i)));
    %Ra=extraterrestrial radiation [MJ m^-2 day^-1] for sub-daily scale
    %Ra(i)=(12*60/pi)*Gsc*dr(i)*((omega2(i)-omega1(i))*sin(latrad)*sin(de(i))+...
     %   cos(latrad)*cos(de(i))*(sin(omega2(i))-sin(omega1(i))))*24;
    %N=daylight hours
    N(i)=24*omegas(i)/3.14;
    %n=actual duration of sunshine, when there are no clouds n=N otherwise n<N
    n(i)=N(i);
    %Rs=solar radiation [MJ m^-2 day^-1]
    Rs(i)=(0.25+0.50*(n(i)/N(i)))*Ra(i);
    %alpha=albedo or canopy reflection coeff. It's 0.23 for the hypothetical grass reference crop
    alpha=0.23;
    %Rns=net solar radiation [MJ m^-2 day^-1]
    Rns(i)=Rs(i)*(1-alpha);
    %Rs0= clear sky solar radiation [MJ m^-2 day^-1]
    Rs0(i)=Ra(i)*(0.75+z*2*(10^(-5)));
    %Rnl=Net longwave radiation [MJ m^-2 day^-1]
    Rnl(i)=4.903*(10^(-9))*(((Tmax(i)+273.16)^(4)+(Tmin(i)+273.16)^(4))/2)*...
        (0.34-0.14*sqrt(ea(i)))*...
        (1.35*(Rs(i)/Rs0(i))-0.35);
    %Rn=Net radiation [MJ m^-2 day^-1]
    Rn(i)=Rns(i)-Rnl(i);
    %G=soil heat flux [MJ m^(-2)]. NEGLECTABLE FOR DAILY SCALES
    G=0;
    %G for hourly scale
%       if omega(i)<-omegas(i) | omega(i)>omegas(i)
%           G=0.5*Rn(i);
%       else
%           G=0.1*Rn(i);
%       end
    
    PETPM(i)=(0.408*Delta(i)*(Rn(i)-G)+gam(i)*(900/(T(i)+273))*u2(i)*(es(i)-ea(i)))/...
        (Delta(i)+gam(i)*(1+0.34*u2(i)));
    % PETPT=Potential evapotranspiration [mm] with Priestley-Taylor
    PETPT(i)=Delta(i)*1.26*(Rn(i)-G)/(lam*(Delta(i)+gam(i)));
    %Evaporation from open water using Shuttleworth (1993).  I multiply by
    %1000 to change the units from m to mm
    EV(i)=1000*(Rn(i)*Delta(i)+6.43*gam(i)*(1+0.536*u2(i))*(es(i)-ea(i)))/...
        (rhow*lam*(Delta(i)+gam(i)));
    
end





end
