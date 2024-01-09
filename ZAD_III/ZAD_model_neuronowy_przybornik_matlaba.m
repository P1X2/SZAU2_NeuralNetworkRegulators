clear all 

dane_wer = readmatrix('dane_wer.txt');
u_wer = dane_wer(:, 1);
y_wer = dane_wer(:, 2);

dane_ucz = readmatrix('dane_ucz.txt');
u_ucz = dane_ucz(:, 1);
y_ucz = dane_ucz(:, 2);

ukm4 = u_ucz(6:length(u_ucz));
ukm3 = u_ucz(7:length(u_ucz));
ykm1 = y_ucz(9:length(u_ucz));
ykm2 = y_ucz(8:length(u_ucz));

ukm4(3996:4000) = 0;
ukm3(3995:4000) = 0;
ykm1(3993:4000) = 0;
ykm2(3994:4000) = 0;
traning_input = {ukm4; ukm3; ykm1; ykm2};

%alguczenia='traingd';%alg. najszybszego spadku
alguczenia='trainlm';%alg. Levenberga-Marquardta
%alguczenia='trainbfg';%alg. zmiennej metryki BFGS
%alguczenia='traincgf';%alg. gradientów sprzężonych Fletchera-Reevesa
%alguczenia='traincgp';%alg. gradientów sprzężonych Poljaka-Polaka-Ribiery


K=6;%liczba neuronów ukrytych
net=feedforwardnet(K,alguczenia);
net.numinputs = 4;
net.configure(net, traning_input);


%sn.performFcn ='mae';%suma modułów błędów/liczba próbek
%sn.performFcn ='mse';%suma kwadratów błędów/liczba próbek
net.performFcn ='sse';
net.trainParam.show = 10;
net.trainParam.showCommandLine = 1;
net.trainParam.epochs = 500;
net.trainParam.goal = 0.0001;
net.trainParam.showWindow = 1;


%dane: tylko zbiór uczący
net.divideFcn = 'divideind';
net.divideParam.trainInd = 1:length(y_ucz);
net.divideParam.valInd = [];
net.divideParam.testInd = [];

net.input.processFcns = { };
net.output.processFcns= { };

[net,uczenie]=train(net,traning_input,y_ucz);

ymod_ucz=sim(net,x_ucz);

Eucz=(y_ucz-ymod_ucz)*(y_ucz-ymod_ucz)';






semilogy(uczenie.perf,'b');
xlabel('Iteracje uczące');
ylabel('Eucz');

figure;
plot(x_ucz,y_ucz,'.b','MarkerSize',14);
hold on;
plot(x_ucz,ymod_ucz,'or');
xlabel('x');
ylabel('y');
legend('Dane','Model');
title(sprintf('Dane uczące, Eucz=%e',Eucz))

figure;
plot(y_ucz,ymod_ucz,'.b','MarkerSize',14);
xlabel('Dane uczące');
ylabel('Model');
title(sprintf('Eucz=%e',Eucz));


% load dane_wer;
% 
% ymod_wer=sim(net,x_wer);
% 
% Ewer=(y_wer-ymod_wer)*(y_wer-ymod_wer)';
% 
% figure;
% plot(x_wer,y_wer,'.b','MarkerSize',14);
% hold on;
% plot(x_wer,ymod_wer,'or');
% xlabel('x');
% ylabel('y');
% legend('Dane','Model');
% title(sprintf('Dane weryfikujące, Ewer=%e',Ewer))
% 
% figure;
% plot(y_wer,ymod_wer,'.b','MarkerSize',14);
% xlabel('Dane weryfikujące');
% ylabel('Model');
% title(sprintf('Ewer=%e',Ewer));



