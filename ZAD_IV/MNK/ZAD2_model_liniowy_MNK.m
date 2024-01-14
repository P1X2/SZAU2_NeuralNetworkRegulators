clear all 

dane_wer = readmatrix('dane_wer.txt');
u_wer = dane_wer(:, 1);
y_wer = dane_wer(:, 2);

dane_ucz = readmatrix('dane_ucz.txt');
u_ucz = dane_ucz(:, 1);
y_ucz = dane_ucz(:, 2);

%% 
dane_ucz = true; % wybór danych
wykresy = true; % włączanie/wyłączanie wykresów
ARX = false;     % wybór typu modelu 

%%
if dane_ucz
    u = u_ucz;
    y = y_ucz;
else
    u = u_wer;
    y = y_wer;
end

%% MNK

steps = length(dane_wer);

Y = y_ucz;
M = zeros(steps, 4);

for i=10:steps
    M(i, :) = [u_ucz(i-3) u_ucz(i-4) y_ucz(i-1) y_ucz(i-2)];
end

w = M\Y;

%%
y_mod = zeros(1, steps);
e = zeros(1, steps);

if ~ARX
    for k=10:steps
        y_mod(k) = w(1) * u(k-3) + w(2) * u(k-4) + w(3) * y_mod(k-1) + w(4) * y_mod(k-2);
        e(k) = y_mod(k) - y(k);
    end

else
    for k=10:steps
        y_mod(k) = w(1) * u(k-3) + w(2) * u(k-4) + w(3) * y(k-1) + w(4) * y(k-2);
        e(k) = y_mod(k) - y(k);
    end

end

Error = sum(e.^2);
disp(Error);

if wykresy
    fig1 = figure;
    hold on 
    plot(y, 'DisplayName', 'y');
    plot(y_mod, '--', "DisplayName", 'y_m_o_d');
    xlabel('k');
    ylabel('y, y_m_o_d')
    error = strcat('Error = ', int2str(Error));
    if ARX
        title("Model ARX" + newline + error)
    else
        title("Model OE" + newline + error)
    end
    legend('Location','southeast')
end










