function Rmse=rmse(obs,mod)
 tmp1=0;
 
 %Ens=sum(diag((obs-mod)*(obs-mod)'))/(sum(obs-mean(obs)))
 for i = 1:length(obs)
 tmp1=tmp1+(obs(i)-mod(i))^2;
 end
 Rmse=sqrt(tmp1/length(obs));