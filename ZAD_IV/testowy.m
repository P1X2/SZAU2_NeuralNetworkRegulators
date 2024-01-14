
% load("OE_BFGS\Najlepsze_wagi_OE_BFGS\N6.mat")
load("MNK\MNK_współczynniki.mat")

% rand('seed', 23);
rand('seed', 70);

alfa1 =  -1.262719;
alfa2 = 0.329193;
beta1 = 0.039291;
beta2 = 0.027184;

    
u = 2 * rand(1, 100) - 1;
U(1:80) = 0;
i = 1;

steps = 4000;
x1 = zeros(1,steps);
x2 = zeros(1,steps);
y = zeros(1,steps);
y_mod = zeros(1,steps);

   
% Symulacja obiektu dla losowych skoków sterowań
for k=10:steps
    if mod(k, 80) == 0 && k ~= 4000 
        U(k+1:k+80) = u(i); 
        i = i+1;
    end

    x1(k) = -alfa1 * x1(k-1) + x2(k-1) + beta1 * g1(U(k-3));
    x2(k) = -alfa2 * x1(k-1) + beta2 * g1(U(k-3));
    y(k) = g2(x1(k));



%     q = [U(k - 3); U(k -4); y(k-1); y(k-2)];
%     y_mod(k) = w20 + w2 * tanh(w10 + w1 * q);



    y_mod(k) = w(1) * U(k-3) + w(2) * U(k-4) + w(3) * y_mod(k-1) + w(4) * y_mod(k-2);
    e(k) = y_mod(k) - y(k);

end

hold on
plot(y, 'DisplayName','y')
plot(y_mod, '--', 'DisplayName', 'y_{mod}')
