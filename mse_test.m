function mse_calc = mse_test(x, net, inputs, targets, nn)
% 'x' contains the weights and biases vector
% in row vector form as passed to it by the
% genetic algorithm. This must be transposed
% when being set as the weights and biases
% vector for the network.
% To set the weights and biases vector to the
% one given as input
%net = setwb(net, x');
net.LW{4,3} = vec2mat(x(1:1*nn(1)),nn(1));
net.LW{3,2} = vec2mat(x(1*nn(1)+1:1*nn(1)+nn(1)*nn(2)),nn(2));
net.LW{2,1} = vec2mat(x(nn(1)*nn(2)+1*nn(1)+1:1*nn(1)+nn(1)*nn(2)+nn(2)*nn(3)),nn(3));

%setwb(net,x);

% To evaluate the ouputs based on the given
% weights and biases vector
y = sim(net, inputs);
%figure(2);
%plot(inputs, targets, 'g'),hold on,plot(inputs, y, 'r'),hold off;pause(10e-10);

% Calculating the mean squared error
%size(y),size(targets),
%r = sum((y-targets).^2)/length(y);
%size(r),

mse_calc = perform(net, targets, y);
%figure(2);
%plot(inputs, targets, 'g'),hold on,plot(inputs, y, 'r'),hold off;pause(10e-10)

end