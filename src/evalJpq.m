% Devuelve las matrices Jp y Jq, usadas para el cáclulo de puntos
% singulares
% Se tiene que cumplir J = Jq^-1*Jp
%
% Si Jp = 0, tenemos singularidad directa (se pierde rigidez)
% Si Jq = 0, tenemos sigunlaridad inversa (se pierde movilidad)
function [Jp, Jq] = evalJpq (miRobot)
    
    % Variables locales de la función
    miq = (miRobot.q)';
    miL = miRobot.L;
    mipsi = [miRobot.psi(1:2); -miRobot.psi(3)];
    beta = miq + mipsi;
    miang = wrapTo2Pi(miRobot.Plataforma.angLcdg' - pi);
    
    % Longitudes proyectadas
    a = miL(:,1) .* [cos(miq)   sin(miq)   [0;0;0]];
    b = miL(:,2) .* [cos(beta)  sin(beta)  [0;0;0]];
    e = miL(:,3) .* [cos(miang) sin(miang) [0;0;0]];
    
    %% Jp
    Jp = b;
    aux = cross(e, b, 2);
    Jp(:,3) = aux(:,3);
    
    %% Jq
    aux = cross(a, b, 2);
    Jq = diag(aux(:,3));    
    
end