load('./trainData/allPoints.mat', 'puntos')
 
% Distancia de cada punto a uno dado-
ptoRef = [0.25 0.22];
for i = 1:length(puntos)
 puntos(i,10) = norm(puntos(i,1:2) - ptoRef);
end
dist = 10;
f = puntos(:,10) < dist;
puntos = puntos(f,:);

% Elegimos la configuración
conf = [1 -1 1];
g = puntos(:,4) == conf(1) & puntos(:,5) == conf(2) & puntos(:,6) == conf(3);
puntos = puntos(g,:);
  
net = feedforwardnet(25);

net.trainParam.showCommandLine = true;
net.trainParam.max_fail = 10;
net.trainParam.lr = 1e-6;
net.trainParam.mc = 0.7;

net.divideParam.trainRatio = 0.7;
net.divideParam.valRatio= 0.15;
net.divideParam.testRatio= 0.15;

% Datos de entrada
puntosq = puntos(:,7:9);
puntosP = puntos(:,1:3);

x = puntosq';
y = puntosP';

% Transformarmos las posiciones para homogeneizar errores
initialIntX = [min(y(1,:)) max(y(1,:))];
initialIntY = [min(y(2,:)) max(y(2,:))];
initialIntPhi = [min(y(3,:)) max(y(3,:))];
newInt = [-1 1];
newIntPhi = newInt/6;
net = configure(net, x, y);
net.name = 'Robot 3RRR';
[net, tr] = train(net, x, y);

save(['nets/', char(datetime("today")), '-net_', num2str(conf, '%d'),'.mat'], 'net');