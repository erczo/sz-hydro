function R2=r2(obs,mod)
 tmp1=0;
 tmp2=0;
 tmp3=0;
 %Ens=sum(diag((obs-mod)*(obs-mod)'))/(sum(obs-mean(obs)))
 for i = 1:length(obs)
 tmp1=tmp1+(obs(i)-mean(obs))*(mod(i)-mean(mod));
 tmp2=tmp2+(obs(i)-mean(obs))^2;
 tmp3=tmp3+(mod(i)-mean(mod))^2;
 end
 R2=tmp1/(sqrt(tmp2)*sqrt(tmp3));