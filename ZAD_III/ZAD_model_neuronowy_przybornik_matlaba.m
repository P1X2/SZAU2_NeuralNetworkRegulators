clear all
%% !!!! najpierw należy przeprowadzić uczenie przy pomocy danych uczacych -> dane = "ucz"!!!!
dane = 'ucz'; % wer, ucz

if dane == "ucz"
    load("dane_ucz.txt")
else
    load("dane_wer.txt")
end

% deklaracja sieci, jej zmiennych i ustawień
input_delay = [3 4];
output_delay = [1 2];
neuron_number = 6;
net = narxnet(input_delay, output_delay, neuron_number);



net.trainFcn = 'trainlm';
narx_net.divideFcn = '';
net.trainParam.epochs = 400;
net.layers{1}.transferFcn = 'tansig';


if strcmp(dane, "ucz")
    X = dane_ucz(:, 1);
    Y = dane_ucz(:, 2);
    X = tonndata(X,false,false);
    Y = tonndata(Y,false,false);
    [Xs,Xi,Ai,Ts] = preparets(net,X,{},Y);
    net = train(net, Xs, Ts, Xi, Ai);  
    Y_pred = sim(net, Xs, Xi, Ai);

    Y_pred = cell2mat(Y_pred);
    Y = cell2mat(Y);
    Y = Y(1:3996);

    Error = sum((Y - Y_pred).^2);

    hold on;
    error = strcat('Error = ', int2str(Error));
    plot(Y_pred, '-', 'DisplayName',"y_m_o_d")
    plot(Y,'DisplayName',"y")
    title("Zbiór uczący" + newline + error)
    xlabel("k")
    ylabel("y")
    legend(Location="northeast")
    legend('show');
end


if strcmp(dane, "wer")
    X = dane_wer(:, 1);
    Y = dane_wer(:, 2);
    X = tonndata(X,false,false);
    Y = tonndata(Y,false,false);
    [Xs,Xi,Ai,Ts] = preparets(net,X,{},Y);
    Y_pred = sim(net, Xs, Xi, Ai);
    
    Y_pred = cell2mat(Y_pred);
    Y = cell2mat(Y);
    Y = Y(1:3996);    

    Error = sum((Y - Y_pred).^2);

    hold on;
    error = strcat('Error = ', int2str(Error));
    plot(cell2mat(Y_pred), '-', 'DisplayName',"y_m_o_d")
    plot(cell2mat(Y),'DisplayName',"y")
    title("Zbiór weryfikujący" + newline + error)
    xlabel("k")
    ylabel("y")
    legend(Location="northeast")
    legend('show');
end
