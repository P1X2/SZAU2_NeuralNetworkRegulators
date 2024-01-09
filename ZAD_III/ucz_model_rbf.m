%uczenie statycznego modelu neuronowego o jednym wejœciu i jednym wyjœciu: 
%sieæ RBF o jednej warstwie ukrytej, automatyczny dobór liczby neuronów
%ukrytych dla b³edu dopuszczalnego (edopuszcz)
clear all;

load dane_ucz;

edopuszcz=0.000000001;

%net = newrb(P,T,goal,spread,MN,DF)
%performFcn
sn=newrb(x_ucz,y_ucz,edopuszcz,1,30,1);%,'performFcn','sse');

ymod_ucz=sn(x_ucz);

Eucz=(ymod_ucz-y_ucz)*(ymod_ucz-y_ucz)';

figure;
plot(x_ucz,y_ucz,'.b','MarkerSize',14);
hold on;
plot(x_ucz,ymod_ucz,'or');
xlabel('x');
ylabel('y');
title(sprintf('Dane ucz¹ce, Eucz=%e',Eucz))
legend('Dane ucz¹ce','Wyjœcie modelu')

figure;
plot(y_ucz,ymod_ucz,'.b','MarkerSize',14);
xlabel('Dane ucz¹ce');
ylabel('Model');
title(sprintf('Eucz=%e',Eucz))

load dane_wer;

ymod_wer=sn(x_wer);

Ewer=(ymod_wer-y_wer)*(ymod_wer-y_wer)';

figure;
plot(x_wer,y_wer,'.b','MarkerSize',14);
hold on;
plot(x_wer,ymod_wer,'or');
xlabel('x');
ylabel('y');
title(sprintf('Dane weryfikuj¹ce, Ewer=%e',Ewer))
legend('Dane weryfikuj¹ce','Wyjœcie modelu')

figure;
plot(y_wer,ymod_wer,'.b','MarkerSize',14);
xlabel('Dane weryfikuj¹ce');
ylabel('Model');
title(sprintf('Ewer=%e',Ewer))





