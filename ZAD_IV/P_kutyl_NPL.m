
clear;
steps_sym=1000;
Upp = 0;
Ypp =0;
u(1:steps_sym)=Upp;
y(1:steps_sym)=Ypp;
e(1:steps_sym) = 0;
q = zeros(steps_sym,4);
q4= zeros(steps_sym,4);
q3= zeros(steps_sym,4);
q2 = zeros(steps_sym,4);
q1 = zeros(steps_sym,4);
%% Przykładowe wartości zadanej yzad
yzad(1:260)=-10;
yzad(261:451)= -1;
yzad(452:700) = -7;
yzad(701:1000) = 0.5;
alfa1 =  -1.262719;
alfa2 = 0.329193;
beta1 = 0.039291;
beta2 = 0.027184;
x1(1:steps_sym) = 0;
x2(1:steps_sym) = 0;
load('OE_BFGS\Najlepsze_wagi_OE_BFGS\N6.mat')
delta = 1e-3;
d(1:steps_sym) = 0;
teta = 3;
y_ob(1:steps_sym) = 0;
D = 100;
N = 300;
Nu=60;
%% Pętla regulatora
lambda = 1000;
Umax = 1;
Umin = -1;

for k=7:steps_sym

    x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(k-3));
    x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(k-3));
    y_ob(k) = g2(x1(k));


    q(k,:) =  [u(k-teta) u(k-teta-1) y_ob(k-1) y_ob(k-2)];

    q1(k,:) = [u(k-teta)+delta u(k-teta-1) y_ob(k-1) y_ob(k-2)];
    q2(k,:) = [u(k - teta) u(k-teta-1)+delta y_ob(k-1) y_ob(k-2)];
    q3(k,:) = [u(k-teta) u(k-teta-1) y_ob(k-1)+delta y_ob(k-2)];
    q4(k,:) = [u(k-teta) u(k-teta-1) y_ob(k-1) y_ob(k-2)+delta];

%     wyjscie ob
    y(k) = w20 + w2*tanh(w10 + w1*q(k,:)'); 

    b3 =  (w20 + w2*tanh(w10 + w1*q1(k,:)') - y(k))/delta;
    b4 =  (w20 + w2*tanh(w10 + w1*q2(k,:)') - y(k))/delta;
    b = [0,0,b3,b4];

    a1 =  - (w20 + w2*tanh(w10 + w1*q3(k,:)') - y(k))/delta;
    a2 =  - (w20 + w2*tanh(w10 + w1*q4(k,:)') - y(k))/delta;
    a = [a1,a2];

    d(k) = y_ob(k) - y(k);


    Y_swobodne(1:N) = 0 ;
    for i=1:N
        if i>=3
            q_pred = [u(k+min(-1,-teta+i)) u(k+min(-1,-teta+i-1)) Y_swobodne(i-1) Y_swobodne(i-2)];
        elseif i==2
            q_pred = [u(k+min(-1,-teta+i)) u(k+min(-1,-teta+i-1)) Y_swobodne(i-1) y_ob(k)];
        else
            q_pred = [u(k+min(-1,-teta+i)) u(k+min(-1,-teta+i-1)) y_ob(k-1+i) y_ob(k-2+i)];
        end
        Y_swobodne(i) = w20 + w2*tanh(w10 + w1*q_pred')+ d(k);
    end

    s(1:N) = 0;
    for j=1:N
        b_czlon = 0;
        a_czlon = 0;
        for i=1:min(j,length(b))
            b_czlon = b_czlon + b(i);
        end
        for p = 1:min(j-1,length(a))
             a_czlon = a_czlon + (a(p) * s(j-p));
        end
        s(j) = b_czlon - a_czlon;
    end

    M = zeros(N,Nu);
    % Macierz M
        for i = 1:N
            for j = 1:Nu
                if (i-j+1) > 0
                    M(i,j) = s(i-j+1);
                end
            end
        end

        Alpha = eye(Nu, Nu) * lambda;
        Yzadk = yzad(k) * ones(N, 1);
        K = inv(M' * M + Alpha) * M';
        dU = K * (Yzadk - Y_swobodne');
        u(k) = dU(1) + u(k-1);
        u(k) = max(min(u(k), Umax), Umin);
        e(k) = (yzad(k) - y_ob(k));
end

%% Przygotowanie wykresów i wizualizacja
t = linspace(1,steps_sym,steps_sym);
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