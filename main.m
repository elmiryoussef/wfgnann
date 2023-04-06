TEXT=readtable('BECHAR.01.02.2005.05.03.2023.1.0.0.en.utf8.00000000.xls','Range','A:H');

M = [];T = [];
for i=round(size(TEXT,1)*0.8):size(TEXT,1)
    va = [];vt = [];
    %va = [str2num(TEXT{i,1}(1:2));str2num(TEXT{i,1}(4:5));str2num(TEXT{i,1}(7:10));str2num(TEXT{i,1}(12:13))];
    va = [str2num(TEXT{i,1}{1}(1:2));str2num(TEXT{i,1}{1}(4:5));str2num(TEXT{i,1}{1}(12:13))];
    M = [M va];
    vt = [TEXT{i,2};TEXT{i,3};TEXT{i,5};TEXT{i,7}];
    %vt = [TEXT{i,2}];
    T = [T vt];
    display(i);
end,
MF=fliplr(M);
TF=fliplr(T);

if isnan(TF(1,1)) TF(1,1)=TF(1,2);end
if isnan(TF(2,1)) TF(2,1)=TF(2,2);end
if isnan(TF(3,1)) TF(3,1)=TF(3,2);end
if isnan(TF(4,1)) TF(4,1)=TF(4,2);end
for i=2:size(TF,2)
    j=find(isnan(TF(:,i)));
    TF(j,i)=TF(j,i-1);
end

% muX = mean(cat(2,MF),2);
% sigmaX = std(cat(2,MF),0,2);

muT = mean(cat(2,TF),2);
sigmaT = std(cat(2,TF),0,2);

for n = 1:size(MF,2)
%     MF(:,n) = (MF(:,n) - muX) ./ sigmaX;
    %TF(:,n) = (TF(:,n) - muT) ./ sigmaT;
end

nn = [10 10 10];
net = feedforwardnet(nn);
n=0;mx = [];
%while n<1
netT = configure(net,MF,TF(1,:));
% netP = configure(netP,MF,TF(2,:));
% netH = configure(netH,MF,TF(3,:));
% netW = configure(netW,MF,TF(4,:));
% Y= sim(net, MF);
% X = perform(net,TF,Y);
% x = net.LW{4,3}';
% if X <= 5e2 mx = [mx;x(:)'];n = n+1,end
% mx = [mx;x(:)'];
%end
%net = configure(net,M,T);
%net.divideFcn = '';
netT.trainParam.showCommandLine = true;
netT.trainParam.show = 100;
netT.trainParam.showWindow = false;
netT.trainParam.max_fail = 1000;
netT.trainParam.epochs = 10000;
netP = netT; netH = netP; netW = netH;
netT1 = netT; netP1 = netP; netH1 = netH; netW1 = netW;

hT = @(x) mse_test(x, netT, MF, TF(1,:), nn);
hP = @(x) mse_test(x, netP, MF, TF(2,:), nn);
hH = @(x) mse_test(x, netH, MF, TF(3,:), nn);
hW = @(x) mse_test(x, netW, MF, TF(4,:), nn);
ga_opts = gaoptimset('InitialPopulation',mx,'display','iter','Generations',10,'TolFun',0,'UseParallel', false);
tic;[x_ga_optT, err_ga] = ga(hT, 1*nn(1) + nn(1)*nn(2) + nn(2)*nn(3), ga_opts);Tga=toc,
tic;[x_ga_optP, err_ga] = ga(hP, 1*nn(1) + nn(1)*nn(2) + nn(2)*nn(3), ga_opts);Tga=toc,
tic;[x_ga_optH, err_ga] = ga(hH, 1*nn(1) + nn(1)*nn(2) + nn(2)*nn(3), ga_opts);Tga=toc,
tic;[x_ga_optW, err_ga] = ga(hW, 1*nn(1) + nn(1)*nn(2) + nn(2)*nn(3), ga_opts);Tga=toc,
%tic;[x_ga_opt, err_ga] = ga(h, nn, ga_opts);Tga=toc,

netT1.LW{4,3} = vec2mat(x_ga_optT(1:1*nn(1)),nn(1));
netT1.LW{3,2} = vec2mat(x_ga_optT(1*nn(1)+1:1*nn(1)+nn(1)*nn(2)),nn(2));
netT1.LW{2,1} = vec2mat(x_ga_optT(nn(1)*nn(2)+1*nn(1)+1:1*nn(1)+nn(1)*nn(2)+nn(2)*nn(3)),nn(3));

netP1.LW{4,3} = vec2mat(x_ga_optP(1:1*nn(1)),nn(1));
netP1.LW{3,2} = vec2mat(x_ga_optP(1*nn(1)+1:1*nn(1)+nn(1)*nn(2)),nn(2));
netP1.LW{2,1} = vec2mat(x_ga_optP(nn(1)*nn(2)+1*nn(1)+1:1*nn(1)+nn(1)*nn(2)+nn(2)*nn(3)),nn(3));

netH1.LW{4,3} = vec2mat(x_ga_optH(1:1*nn(1)),nn(1));
netH1.LW{3,2} = vec2mat(x_ga_optH(1*nn(1)+1:1*nn(1)+nn(1)*nn(2)),nn(2));
netH1.LW{2,1} = vec2mat(x_ga_optH(nn(1)*nn(2)+1*nn(1)+1:1*nn(1)+nn(1)*nn(2)+nn(2)*nn(3)),nn(3));

netW1.LW{4,3} = vec2mat(x_ga_optW(1:1*nn(1)),nn(1));
netW1.LW{3,2} = vec2mat(x_ga_optW(1*nn(1)+1:1*nn(1)+nn(1)*nn(2)),nn(2));
netW1.LW{2,1} = vec2mat(x_ga_optW(nn(1)*nn(2)+1*nn(1)+1:1*nn(1)+nn(1)*nn(2)+nn(2)*nn(3)),nn(3));

tic;[netT1, tr1] = train(netT1,MF,TF(1,:),'useGPU','yes');Tnn1=toc,
%net.LW{4,3} = vec2mat(mean(mx,1),nn(1));
tic;[netT, tr] = train(netT,MF,TF(1,:),'useGPU','yes');Tnn=toc,

tic;[netP1, tr1] = train(netP1,MF,TF(2,:),'useGPU','yes');Tnn1=toc,
%net.LW{4,3} = vec2mat(mean(mx,1),nn(1));
tic;[netP, tr] = train(netP,MF,TF(2,:),'useGPU','yes');Tnn=toc,

tic;[netH1, tr1] = train(netH1,MF,TF(3,:),'useGPU','yes');Tnn1=toc,
%net.LW{4,3} = vec2mat(mean(mx,1),nn(1));
tic;[netH, tr] = train(netH,MF,TF(3,:),'useGPU','yes');Tnn=toc,

tic;[netW1, tr1] = train(netW1,MF,TF(4,:),'useGPU','yes');Tnn1=toc,
%net.LW{4,3} = vec2mat(mean(mx,1),nn(1));
tic;[netW, tr] = train(netW,MF,TF(4,:),'useGPU','yes');Tnn=toc,

Mt = [];Tt = [];
for i=2:round(size(TEXT,1)*0.2)
    vat = [];vtt = [];
    %vat = [str2num(TEXT{i,1}(1:2));str2num(TEXT{i,1}(4:5));str2num(TEXT{i,1}(7:10));str2num(TEXT{i,1}(12:13))];
    vat = [str2num(TEXT{i,1}{1}(1:2));str2num(TEXT{i,1}{1}(4:5));str2num(TEXT{i,1}{1}(12:13))];
    Mt = [Mt vat];
    vtt = [TEXT{i,2};TEXT{i,3};TEXT{i,5};TEXT{i,7}];
    %vtt = [TEXT{i,2}];
    Tt = [Tt vtt];
end,
if isnan(Tt(1,1)) Tt(1,1)=Tt(1,2);end
if isnan(Tt(2,1)) Tt(1,1)=Tt(2,2);end
if isnan(Tt(3,1)) Tt(1,1)=Tt(3,2);end
if isnan(Tt(4,1)) Tt(1,1)=Tt(4,2);end
for i=2:size(Tt,2)
    j=find(isnan(Tt(:,i)));
    Tt(j,i)=Tt(j,i-1);
end
%rt(find(isnan(Tt)))=0;
%Tt(find(isnan(Tt)))=0;

Mt=fliplr(Mt);
Tt=fliplr(Tt);

% muX = mean(cat(2,Mt),2);
% sigmaX = std(cat(2,Mt),0,2);

muT = mean(cat(2,Tt),2);
sigmaT = std(cat(2,Tt),0,2);

for n = 1:size(Mt,2)
%     Mt(:,n) = (Mt(:,n) - muX) ./ sigmaX;
    %Tt(:,n) = (Tt(:,n) - muT) ./ sigmaT;
end

rt1(1,:) = sim(netT1,Mt);
rt1(2,:) = sim(netP1,Mt);
rt1(3,:) = sim(netH1,Mt);
rt1(4,:) = sim(netW1,Mt);

rmse1 = sqrt(mean((Tt-rt1).^2,2))',
rr1 = abs(Tt-rt1);
erm1 = mean(rr1');
ermx1 = max(rr1');
ermn1 = min(rr1');

rt(1,:) = sim(netT,Mt);
rt(2,:) = sim(netP,Mt);
rt(3,:) = sim(netH,Mt);
rt(4,:) = sim(netW,Mt);

rmse = sqrt(mean((Tt-rt).^2,2))',
rr = abs(Tt-rt);
erm = mean(rr');
ermx = max(rr');
ermn = min(rr');

subplot(2,2,1),plot(1:round(size(Tt,2)), Tt(1,1:round(size(Tt,2))),'--r'),hold on, plot(1:round(size(rt,2)), rt(1,1:round(size(rt,2))),'--b'), plot(1:round(size(rt1,2)), rt1(1,1:round(size(rt1,2))),'--g'), hold off,
xlabel('Time'), ylabel('Air temperature (Degrees celsius)'), legend('Measured data','ANN data','GA+ANN data'), title('Air temperature curves.')
subplot(2,2,2),plot(1:round(size(Tt,2)), Tt(2,1:round(size(Tt,2))),'--r'),hold on, plot(1:round(size(rt,2)), rt(2,1:round(size(rt,2))),'--b'), plot(1:round(size(rt1,2)), rt1(2,1:round(size(rt1,2))),'--g'),hold off,
xlabel('Time'), ylabel('Atmospheric pressure (Millimeters of mercury)'), legend('Measured data','ANN data','GA+ANN data'), title('Atmospheric pressure curves.')
subplot(2,2,3),plot(1:round(size(Tt,2)), Tt(3,1:round(size(Tt,2))),'--r'),hold on, plot(1:round(size(rt,2)), rt(3,1:round(size(rt,2))),'--b'), plot(1:round(size(rt1,2)), rt1(3,1:round(size(rt1,2))),'--g'), hold off,
xlabel('Time'), ylabel('Relative humidity (%)'), legend('Measured data','ANN data','GA+ANN data'), title('Relative humidity curves.')
subplot(2,2,4),plot(1:round(size(Tt,2)), Tt(4,1:round(size(Tt,2))),'--r'),hold on, plot(1:round(size(rt,2)), rt(4,1:round(size(rt,2))),'--b'), plot(1:round(size(rt1,2)), rt1(4,1:round(size(rt1,2))),'--g'), hold off,
xlabel('Time'), ylabel('Mean wind speed (Meters per second)'), legend('Measured data','ANN data','GA+ANN data'), title('Mean wind speed curves.')