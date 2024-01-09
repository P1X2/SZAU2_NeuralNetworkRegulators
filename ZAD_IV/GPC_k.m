
clear all
% wagi najlepszego modelu uczonegi w trybie rekurencyjnym
% load('OE_BFGS\Najlepsze_wagi_OE_BFGS\N6.mat')

%% Wybór regulatora; definicja zmiennych modelu
regulator = "PID"; % PID, GPC, NPL, NO
plots = true;
steps_sym = 660;

alfa1 =  -1.262719;
alfa2 = 0.329193;
beta1 = 0.039291;
beta2 = 0.027184;

%% Definicja trajektorii zadanej i wymaganych macierzy

% trajektoria zadana
k_step1 = 13; 
k_step2 = 180;
k_step3 = 260;
k_step4 = 340;
k_step5 = 420;
k_step6 = 500;
k_step7 = 580;

% yzad
yzad(1:k_step1-1) = 0 ; yzad(k_step1:k_step2) = -3;  % yzad z  przedziału  <-11.5, 0.3> bo ograniczenia na u 
yzad(k_step2-1:k_step3) = -5; yzad(k_step3:k_step4-1) = .5;
yzad(k_step4:k_step5-1) = -2; yzad(k_step5:k_step6-1) = -10;
yzad(k_step6:k_step7-1) = -1; yzad(k_step7:steps_sym) = -7;


u(1:k_step1-1) = 0;
x1 = zeros(1,steps_sym);
x2 = zeros(1,steps_sym);
y = zeros(1,steps_sym);
e = zeros(1, steps_sym);

M = zeros(N,Nu);
s(1:N) = 0;

b5 = 

teta = 3;

tested_y_zadk = []
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
lambda =1160;
Umax =1;
Umin = -1;
y_mod(1:kk) = 0;
Alpha = eye(Nu, Nu) * lambda;
K = inv(M' * M + Alpha) * M';
for k=7:steps_sym
    x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(k-3));
    x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(k-3));
    y(k) = g2(x1(k));
    y_mod(k) = b5 * u(k-teta) + b6 * u(k-teta-1) - a1* y_ob(k-1) - a2 * y_ob(k-2);
    Y_swobodne(1:N) = 0;
    d(k) = y_ob(k) - y_mod(k);
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
matlab2tikz ('zad5dodPID.tex' , 'showInfo' , false)