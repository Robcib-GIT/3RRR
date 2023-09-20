% Evaluar los tres métodos cerca de las singularidades.
% Jorge F. García-Samartín
% wwww.gsamartin.es
% 2023-04-13

% Cargamos todo
load('2022-10-24-net_-1 1-1.mat');
robotNuestro;

% Lista de puntos singulares
if ~exist('puntos', 'var')
    [~, ~, ~, ~, puntos] = r.medirParams(0,[-1 1 -1], 'tolP', .5e-4, 'tolQ', .5e-5, 'res', 0.005);
end
singsList = find(puntos(:,7) ~= 0);

% Evaluamos
T = zeros(length(singsList), 3);
A = zeros(length(singsList), 3);
diary 'evalSings.txt'
for p = 1:length(singsList)
    
    Preal = puntos(singsList(p), 1:2);

    % MCD Numérico
    r.P = [0 0 0];
    t1 = datevec(datetime);
    r.q = puntos(singsList(p), 4:6);
    r.MCD;
    t2 = datevec(datetime);
    T(p,1) = norm(t2-t1);
    A(p,1) = norm(Preal - r.P(1:2));
    
    % MCD Neuronal
    t1 = datevec(datetime);
    q = puntos(singsList(p), 4:6)';
    P = sim(net, q);
    t2 = datevec(datetime);
    T(p,2) = norm(t2-t1);
    A(p,2) = norm(Preal - P(1:2));

    % MCD Genético
    r.P = [0 0 0];
    t1 = datevec(datetime);
    r.q = puntos(singsList(p), 4:6);
    r.MCDGen;
    t2 = datevec(datetime);
    T(p,3) = norm(t2-t1);
    A(p,3) = norm(Preal - r.P(1:2));

end
diary off

% Quitamos los valores atípicos
A2 = rmmissing(A);
mean(A2);
A3 = A2(A2(:,1) < 1000,:);