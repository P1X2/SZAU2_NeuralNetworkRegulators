
clear;
load('MNK_współczynniki.mat')
kk=1000;
Upp = 0;
Ypp =0;
u(1:kk)=Upp;
y(1:kk)=Ypp;
e(1:kk) = 0;
%% Przykładowe wartości zadanej yzad
yzad(1:260)=-1.33;
yzad(261:451)= -10;
yzad(452:700) = -7;
yzad(701:1000) = -2;

alfa1 =  -1.262719;
alfa2 = 0.329193;
beta1 = 0.039291;
beta2 = 0.027184;

x1(1:kk) =0;
x2(1:kk) = 0;
delta = 1e-3;
d(1:kk) = 0;
teta = 3;
y_ob(1:kk) = 0;
D = 1000;
N = 160;
Nu=60;
b5 = W(1);
b6 = W(2);
a1 = -W(3);
a2 = -W(4);
b = [0,0,b5,b6];
a = [a1,a2];
M = zeros(N,Nu);
s(1:N) = 0;


for j=1:N
    for i=1:min(j,length(b))
        s(j) = s(j) + b(i);
    end
    for i = 1:min(j-1,length(a))
         s(j) = s(j) - (a(i) * s(j-i));
    end
end


for i = 1:N
    for j = 1:Nu
        if (i-j+1) > 0
            M(i,j) = s(i-j+1);
        end
    end
end
%% Pętla regulatora
lambda =200160;
Umax =1;
Umin = -1;
y_mod(1:kk) = 0;
Alpha = eye(Nu, Nu) * lambda;
K = ((M' * M + lambda * eye(Nu, Nu))^-1) * M';

for k=7:kk

    x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(k-3));
    x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(k-3));
    y_ob(k) = g2(x1(k));

    y_mod(k) = b5 * u(k-teta) + b6 * u(k-teta-1) - a1* y_ob(k-1) - a2 * y_ob(k-2);
    d(k) = y_ob(k) - y_mod(k);

    Y_swobodne(1:N) = 0;
    for i=1:N
        if i>=3
            Y_swobodne(i) = b5 * u(k+min(-1,-teta+i)) + b6 * u(k+min(-1,-teta+i-1)) - a1* Y_swobodne(i-1) - a2 * Y_swobodne(i-2) + d(k);
        elseif i==2
            Y_swobodne(i) = b5 * u(k+min(-1,-teta+i)) + b6 * u(k+min(-1,-teta+i-1)) - a1* Y_swobodne(i-1) - a2 * y_ob(k) + d(k);
        else
            Y_swobodne(i) = b5 * u(k+min(-1,-teta+i)) + b6 * u(k+min(-1,-teta+i-1)) - a1* y_ob(k-i+1) - a2 * y_ob(k-2+i) + d(k);
        end
    end
    Yzadk = yzad(k) * ones(N, 1);
    dU = K * (Yzadk - Y_swobodne');
    u(k) = dU(1) + u(k-1);
    u(k) = max(min(u(k), Umax), Umin);
    e(k) = (yzad(k) -y_ob(k));
end

%% Przygotowanie wykresów i wizualizacja
t = linspace(1,kk,kk);
figure
subplot(2,1,2)
stairs(t,u,'LineWidth',1.5, Color='r');
title('u - sterowanie');
xlabel('k - number próbki');
ylabel("Wartość sterowania")
subplot(2,1,1)
% matlab2tikz ('zad4PID_u.tex' , 'showInfo' , false)
stairs(t,y_ob,'LineWidth',1.5);
hold on
stairs(t,yzad,'LineWidth',1, 'LineStyle','--');
title('Charakterystyki y,y_{zad}');
xlabel('k - number próbki');
ylabel('Wartość')
legend("Wartość na wyjściu y", "Wartość zadana y_{zad}",Location="northwest")
% matlab2tikz ('zad5dodPID.tex' , 'showInfo' , false)