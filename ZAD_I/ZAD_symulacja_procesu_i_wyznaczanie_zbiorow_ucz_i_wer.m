clear all

zadanie = 2; % 1 - punkt pracy; 2 - charakterystyka statyczna; 3 - generowanie zbioru uczącego i weryfikującego


alfa1 =  -1.262719;
alfa2 = 0.329193;
beta1 = 0.039291;
beta2 = 0.027184;

%% sprawdznie pkt pracy
if zadanie == 1
    steps = 100;
    u = zeros(1,steps);
    u(50:steps) = .1;
    x1 = zeros(1,steps);
    x2 = zeros(1,steps);
    y = zeros(1,steps);
    
    % Sprawdzenie pkt pracy
    for k=10:steps
        x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(k-3));
        x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(k-3));
        y(k) = g2(x1(k));
        
    end

    fig1 = figure;
    hold on
    plot(y)
    xlabel('k')
    ylabel('y')
    title('Punkt pracy')
end


%% wyznaczenie charakterystyki statycznej
if zadanie == 2
    steps = 100;
    
    u = -1:0.01:1;
    ystat = zeros(1, length(u));
    
    for i=1:length(u)
    
        x1 = zeros(1,steps);
        x2 = zeros(1,steps);
        y = zeros(1,steps); 
    
        for k=10:steps
            x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(u(i));
            x2(k) = -alfa2 * x1(k-1) + beta2 * g1(u(i));
            y(k) = g2(x1(k));
        end
    
        ystat(i) = y(steps-1);
    
    end
    
    fig2 = figure;
    plot(u, ystat)
    xlabel('u')
    ylabel('y')
    title('Charakterystyka statyczna')

end

%% generowanie danych uczących i weryfikujących

if zadanie == 3
       
    steps = 4000;
    x1 = zeros(1,steps);
    x2 = zeros(1,steps);
    y = zeros(1,steps);

    %% generowanie danych uczących i weryfikujących

    dane_wer = true; % wybór typu danych do generowania 
    szum_pomiarowy = true;

    if dane_wer
        rand('seed', 23);
    else
        rand('seed', 43);
    end
        
    u = 2 * rand(1, 100) - 1;
    U(1:80) = 0;
    i = 1;
       
    % Symulacja obiektu dla losowych skoków sterowań
    for k=10:steps
        if mod(k, 80) == 0 && k ~= 4000 
            U(k+1:k+80) = u(i); 
            i = i+1;
        end

        x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(U(k-3));
        x2(k) = -alfa2 * x1(k-1) + beta2 * g1(U(k-3));
        if szum_pomiarowy
            y(k) = g2(x1(k)) + (0.03 * (2 * rand(1,1) - 1)); % y z szumem pomiarowym 
        else
            y(k) = g2(x1(k));
        end

    end
    
    %% wykresy
    fig1 = figure;
    plot(y)
    xlabel('k')
    ylabel('y')
    if dane_wer
        title('Dane weryfikujące');
    else
        title('Dane uczące')
    end
    

    %% zapis do pliku 
    if dane_wer
        dane = [U', y'];
        fileID = fopen("dane_wer.txt", 'w');
        fprintf(fileID, '%f\t%f\n', dane');
        fclose(fileID);
    else
        dane = [U', y'];
        fileID = fopen("dane_ucz.txt", 'w');
        fprintf(fileID, '%f\t%f\n', dane');
        fclose(fileID);
    end
    
end







