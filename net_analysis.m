function [link_no,downstream_link_no,upstream_link_no1,upstream_link_no2,order,len,...
    magnitude,downstream_contr_area,slope,upstream_contr_area...
    ,watershed_no,watershed_area,wid,new_order]= net_analysis(c,b,net1,tree);

%Make sure the input is ordered properly

%for the Elder use flag1=1, for the Eel river use flag1=0
flag1=1;
% if flag1
% csvread net150.csv;   % <-- moved loading to main code, cab.
% load elder30mtree150.dat
% tree=elder30mtree150;
% else
% csvread Eelnet.csv;
% load Eeltree.dat
% tree=Eeltree; 
% end


%net1=ans;
f= net1(:,9) > 0;
net1=net1(f,:);
clear ans
if flag1
    for i = 1:length(net1(:,1))
        net1(i,3)=tree(length(net1(:,1))-i+1,1);
    end
else
    net1(:,3)=net1(:,4);
end
outlet=find(net1(:,5)==-1);
net(length(net1(:,1)),:)=net1(outlet,:);
i=length(net1(:,1))-1;
j=length(net1(:,1))-1;


%CHECK UNITS IN THE NETWORK FILE PRODUCED BY MAPWINDOW OR ARCGIS
while j>0
     if net(j+1,6) >0
         net(i,:)=net1(find(net1(:,3)==net(j+1,6)),:);
         net(i-1,:)=net1(find(net1(:,3)==net(j+1,7)),:);
         i=i-2;
     end
     j=j-1;
end

for i=1:length(net(:,1))
    if net(i,9)==1
        new_order(i)=1;
    else
        if net(net(:,3)==net(i,6),9)==net(net(:,3)==net(i,7),9)
            new_order(i)=net(net(:,3)==net(i,6),9)+1;
        else
            new_order(i)=max(net(net(:,3)==net(i,6),9),net(net(:,3)==net(i,7),9))+1;
        end
    end    
end

x=net(:,1);
y=net(:,2);
link_no=net(:,3);
downstream_link_no=net(:,5);
upstream_link_no1=net(:,6);
upstream_link_no2=net(:,7);
order=net(:,9);
len=net(:,10); 
magnitude=net(:,11);
downstream_contr_area=net(:,12); %[m^2]
for ii=1:length(len)
    slope(ii)=net(ii,13)/len(ii);
end
upstream_contr_area=net(:,16); %[m^2] 
watershed_no=net(:,17);
watershed_area=net(:,12)-net(:,16);
for ii=1:length(len)
    if net(ii,9)==1
        watershed_area(ii)=net(ii,12);
    end
end
downstream_to_outlet=net(:,18);
upstream_to_outlet=net(:,19);
midstream_to_outlet=net(:,20);
%Channel width using Montgomery and Gran 2001 WRR (Table 1), 
%w=c*A^b, w is in meters and A is in square meters.
%I use an averge contributing area between up- and down-stream: [net(:,12)+net(:,16)]/2
%c=0.01;
%b=0.39;
wid=c*((net(:,12)+net(:,16))/2).^b; %[m]

%save('C:/Users/stefano/Dropbox/Algal model/Algae hydro model/network.mat',...
 %   'link_no','downstream_link_no','upstream_link_no1','upstream_link_no2','order','len',...
 %   'watershed_no','watershed_area');