clear all
% wagi najlepszego modelu uczonegi w trybie rekurencyjnym
% load('OE_BFGS\Najlepsze_wagi_OE_BFGS\N6.mat')

%% Wybór regulatora; definicja zmiennych modelu
regulator = "GPC"; % PID, GPC, NPL, NO
plots = true;
steps_sym = 660;


% Parametry NPL i GPC
N = 30;
Nu = 3;
Lambda = 0.1;



alfa1 = -1.262719;
alfa2 = 0.329193;
beta1 = 0.039291;
beta2 = 0.027184;

%% Definicja trajektorii zadanej i wymaganych macierzy

% trajektoria zadana
k_step1 = 100; 
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp("PID", regulator)) 
    %% Parametry PID-a
    k = 0.1; % to są najlepsze nastawy, ale do sprawka trzeba wrzucić k = += 0.2
    Ti = 12; 
    Td = .1;  
    Tp = 1; 
    
    
    r0 = k*(1+(Tp/(2*Ti))+(Td/Tp));
    r1 = k*((Tp/(2*Ti))-((2*Td)/Tp)-1);
    r2 = (k*Td)/Tp;

    %% Pętla symulacji
    for k=10:steps_sym

        x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(k-3));
        x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(k-3));
        y(k) = g2(x1(k));
        
        e(k) = yzad(k) - y(k);
        u(k) = r2*e(k-2) + r1*e(k-1) + r0*e(k) + u(k-1);
    
        if u(k) <= -1
            u(k) = -1;
        
        elseif u(k) >= 1 
            u(k) = 1;

        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp("GPC", regulator))
    %%  Parametry GPC i NPL znajdują się na początku pliku

    
    % użyto modelu wyznaczonego metodą najmniejszych kwadratów wyznaczonego w
    % zadaniu 2
    load('MNK_współczynniki.mat');
    
    b = [0 0 w(1) w(2)];
    a = [-w(3) -w(4)];
    
    
     
    % wyznaczenie współczynników odpowiedzi skokowej z modelu 
    s = zeros(1,N);
    for i=1:N
        b_sum = 0;
        a_sum = 0;
    
        for ii=1:min(i,4)
            b_sum = b_sum + b(ii);
        end
        
        for jj=1:min(i-1, 2)
            a_sum = a_sum + (a(jj) * s(i-jj));
        end
        
        s(i) = b_sum - a_sum;
    end
     
%     plot(s);
    
    % Macierz M
    M=zeros(N,Nu); 
    for i=1:Nu
        M(i:N,i)=s(1:(N-i+1));
    end
    
    % Współczynnik K
    K = ((M' * M + Lambda * ones(Nu, Nu))^-1) * M';
    
    %% Główna pętla symulacji

    for k=10:steps_sym
    
        % obiekt
        x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(k-3));
        x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(k-3));
        y(k) = g2(x1(k));

        
        % d(k)
        y_mod = b(3) * u(k-3) + b(4) * u(k-4) - a(1) * y(k-1) - a(2) * y(k-2);
        d = y(k) - y_mod;


        % odp. swobodna 
        for i=1:N
            y(k+i) = b(3) * u(min(k-1, k-3+i)) + b(4) * u(min(k-1, k-4+i)) - a(1) * y(k+i-1) - a(2) * y(k+i-2) + d;
        end

        for

        Y_swobodny = y(k+1:k+N);
        Y_zad = ones(N,1) * yzad(k);
        y = y(1:k);

        % wektor du_k
        du = K * (Y_zad - Y_swobodny);


        u(k) = u(k-1) + du(1);
        
        if u(k) > 1
            u(k) = 1;
        elseif u(k) < -1
            u(k) = -1;
        end
        
        e(k) = yzad(k) - y(k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp("NPL", regulator)

    load('OE_BFGS\Najlepsze_wagi_OE_BFGS\N6.mat')
    
    for k=10:steps_sym
    
        % obiekt
        x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(k-3));
        x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(k-3));
        y(k) = g2(x1(k));
    
    
        % d(k)
        q = [u(k - 3); u(k -4); y(k-1); y(k-2)];
        y_mod = w20 + w2 * tanh(w10 + w1 * q);
        d = y(k) - y_mod;

        
        % trajektoria swobodna
        for p=1:N
            q = [u(min(k-1, k-3+i)); u(min(k-1, k-4+i)); y(k+i-1); y(k+i-2)];
            y(k+p) = w20 + w2 * tanh(w10 + w1 * q);
        end
        Y_swobodny = y(k+1:k+N);
        Y_zad = yzad(k) * ones(N, 1);

        
        % odpowiedź skokowa
    
    
    
    
    
    
    end






























end

error_sum = sum(e.^2);

%% Plots
if plots == true

    tit1 = strcat("Error = ", int2str(error_sum));
    fig1=figure;
    subplot(2,1,1);
    hold on
    stairs(y(1:steps_sym), "DisplayName","y")
    stairs(yzad, "DisplayName","y_z_a_d")
    xlabel('k')
    ylabel('y/y_z_a_d')
    legend('Location','southwest')
    title(tit1 + newline + 'Wyjście');
    hold off
    
    
    subplot(2,1,2)
    stairs(u, "DisplayName","u")
    xlabel('k')
    ylabel('u')
    legend('Location','southwest');
    title('Sterowanie');

end


function [] = funreggpc()
global N Nu yzad y u du d k K b a;
    

    
end
