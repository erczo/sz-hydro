function Ens=nashsutcliffe(obs,mod)
 tmp1=0;
 tmp2=0;
 %Ens=sum(diag((obs-mod)*(obs-mod)'))/(sum(obs-mean(obs)))
 for i = 1:length(obs)
 tmp1=tmp1+(obs(i)-mod(i))^2;
 tmp2=tmp2+(obs(i)-mean(obs))^2;
 end
 Ens=1-tmp1/tmp2;