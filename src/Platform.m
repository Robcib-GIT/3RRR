%% Plataforma central del Robot 3RRR 
% 
% TFG de Jorge García Samartín (2020)
%
% Declaración del objeto
%
% CALL plataforma = Platform(T)
%
%   plataforma = objeto con la plataforma
%   T = Longitudes de los lados de la plataforma (vector 1x3)
%
% Ejemplo:
%   plataforma = Platform([2 2 3.5])
%
% La plataforma se construye en ejes locales suponiendo el primer lado
% horizontal y como origen de coordenadas el primer vértice.
% Si posteriormente se define el valor del centro P, se obienen las
% coordendas globales de la plataforma en C
%
%% Esquema de la plataforma triángular
%  
% <<../figuras/esquema.png>>
% 
classdef Platform < handle
% Clase para gestionar el robot
% -------------------------------------------------------------------------
% J.G.S. v1.01 2020-01-30
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Declaración de propiedades visibles accesibles desde fuera de la clase
% -------------------------------------------------------------------------
properties
    T = [0.2 0.2 0.2];      % Longitudes de los lados de la plataforma
    P = [0 0 0];            % Coord. globales impuestas CDG de plataforma
end

% -------------------------------------------------------------------------
% Declaración de propiedades dependientes que pueden asignarse
% -------------------------------------------------------------------------
properties (SetAccess = public, Dependent = true)
    nLados                  % Nº de lados de la plataforma
    equilatero              % 1 si el triángulo de la plataforma es equilátero
    theta = zeros(3,1);     % Angulos interiores 
    coord = zeros(3,2);     % Coordenadas locales vértices de plataforma
    cdg = zeros(1,2);       % Coordenadas del cdg de la plataforma
    L = zeros(3,1);         % Longitud de las líneas del vértice al cdg
    gamma = zeros(3,1);     % Angulos del lado a linea de cdg
    angLcdg = zeros(3,1);   % Ang. global desde horizontal a linea de cdg a B
    B = zeros(1,2);         % Coordenadas globales vértices de plataforma
end
% -------------------------------------------------------------------------
% Declaración de propiedades ocultas accesibles desde fuera de la clase
% -------------------------------------------------------------------------
properties (Hidden, Constant)

end

% -------------------------------------------------------------------------
% Declaración de propiedades ocultas no accesibles desde fuera de la clase
% -------------------------------------------------------------------------
properties (Hidden, SetAccess = private, Dependent = true)
    
end

% -------------------------------------------------------------------------
% Declaración de eventos para detectar cambios
% -------------------------------------------------------------------------
events

end


%% Declaración de métodos
methods

%% Constructor. 
% Método principal de la clase. Se ejecuta cuando se carga un nuevo objeto
%
% Puede llamarse mediante;
%   obj = Platform(T)
%
function obj = Platform(T, P) 
    obj.T = T;
    obj.P = P;
end

%% Salida en pantalla del objeto
% Esta función gestiona la impresión en pantalla del objeto
% function disp(obj)
%     fprintf('Características generales del robot 3RRR:\n');
%     fprintf('%s [%8.4f, %8.4f;\n %13.4f, %8.4f;\n %13.4f, %8.4f]; %s\n', 'O = ', obj.O', ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f;\n %13.4f, %8.4f, %8.4f;\n %13.4f, %8.4f, %8.4f]; %s\n', 'L = ', obj.L', ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'T = ', obj.T, ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 's = ', obj.s, ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'q = ', obj.q, ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'P = ', obj.P, ' ');
%     fprintf('%s [%8.4f, %8.4f;\n %13.4f, %8.4f;\n %13.4f, %8.4f]; %s\n', 'A = ', obj.A', ' ');
%     fprintf('%s [%8.4f, %8.4f;\n %13.4f, %8.4f;\n %13.4f, %8.4f]; %s\n', 'B = ', obj.B', ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'sgn = ', obj.sgn, ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'psi = ', obj.psi, ' ');
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'theta = ', obj.theta, ' ');  % ********
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'gamma = ', obj.gamma, ' ');  % ********
%     fprintf('%s [%8.4f, %8.4f, %8.4f]; %s\n', 'mierror = ', obj.mierror, ' ');
% end

%% Parámetro T. 
function obj = set.T(obj, T)
    % Longitud de los lados de la plataforma
    if isempty(T)
        obj.T = obj.T;
    else
        
        obj.T = T;
%         
%         % Comprobar si es un triángulo
%         if length(T) == 3
%             % Comprobar que la suma de 2 lados es mayor que el tercero
%             if sum(T(1:2)) < T(3) || sum(T(2:3)) < T(1) || sum(T(1:2:3)) < T(2)
%                 error('Tríángulo central mal definido')
%             else      
%                 obj.T = T;
%             end
%         else
%             error('La plataforma no es un triángulo')
%         end
    end
end

function T = get.T(obj)
    T = obj.T;
end


%% Parámetro P. 
function obj = set.P(obj, P)
    % Coordenadas globales impuestas del CDG de la plataforma
    % P(x, y, giro)
    if isempty(P)
        obj.P = obj.P;
    else
        obj.P = P;
    end
end

function P = get.P(obj)
    P = obj.P;
end


%% Propiedades dependientes (solo tienen método get)
% Las propiedades que se definen a continuación son dependientes de los
% datos de entrada y se calculan cada vez que cambiamos uno de ellos. Por
% tanto no se les puede asignar un valor como a los datos de entrada. Ello
% implica que son de solo de lectura, tienen método GET pero no método SET

%% Parámetro nlados
function nLados = get.nLados(obj)
    nLados = length(obj.T);
end

%% Parámetro equilatero
function equilatero = get.equilatero(obj)
    % Devuelve 1 si el triángulo es equilátero y 0 en caso contrario
    
    % Variables locales de la función
    miT = obj.T;
    
    if miT(1) == miT(2) && miT(1) == miT(3)
        equilatero = 1;
    else
        equilatero = 0;
    end
end

%% Parámetro theta
function theta = get.theta(obj)
    % Angulos interiores entre lados de la plataforma
    
    % Variables locales de la función
    miT = obj.T;
    
    if obj.nLados == 3
        if obj.equilatero
            theta = pi/3*[1 1 1];
        else
            theta = [acos(-(miT(2)^2 - miT(1)^2 - miT(3)^2) / 2 / miT(1) / miT(3)), ...
                     acos(-(miT(3)^2 - miT(2)^2 - miT(1)^2) / 2 / miT(2) / miT(1)), ...
                     acos(-(miT(1)^2 - miT(3)^2 - miT(2)^2) / 2 / miT(3) / miT(2))];
        end
    end
end

%% Parámetro c
function coord = get.coord(obj)
    % Coordenadas locales de los vértices de la plataforma.
    % Suponemos el primer punto en 0,0 y el primer lado horizontal
    % matriz con puntos en filas y x,y en columnas
    
    % Variables locales de la función
    miT = obj.T;
    mitheta = obj.theta;
    
    if obj.nLados == 3
        coord = [0, 0; ...
             miT(1), 0; ...
             miT(3) * cos(mitheta(1)), miT(3) * sin(mitheta(1))];
    end
end

%% Parámetro cdg
function cdg = get.cdg(obj)
    % Coordenas x é y del centro de la plataforma
    if obj.nLados == 3
        cdg = mean(obj.coord) ;
    end
end

%% Parámetro L
function L = get.L(obj)
    % Longitud de las líneas del vértice al cdg
    
    if obj.nLados == 3
        if obj.equilatero
            L = obj.T/cos(pi/6)/2;
        else
            % Variables locales de la función
            micoord = obj.coord;
            micdg = obj.cdg;
            
            L = sqrt((micdg(1) - micoord(:,1)).^2 + (micdg(2) - micoord(:,2)).^2);
        end
    end
end

%% Parámetro gamma
function gamma = get.gamma(obj)
    % Angulos interiores entre lado de la plataforma y la recta al cdg
    
    if obj.nLados == 3 
        if obj.equilatero
            gamma = pi/6 * [1 1 1];
        else
            % Variables locales de la función
            miT = obj.T;
            miL = obj.L;

            gamma = [acos(-(miL(2)^2 - miL(1)^2 - miT(1)^2) / 2 / miL(1) / miT(1)), ...
                     acos(-(miL(3)^2 - miL(2)^2 - miT(2)^2) / 2 / miL(2) / miT(2)), ...
                     acos(-(miL(1)^2 - miL(3)^2 - miT(3)^2) / 2 / miL(3) / miT(3))];
        end
    end
end

%% Parámetro angLcdg
function angLcdg = get.angLcdg(obj)
    % Angulo local entre recta horizontal por cdg y L
    if obj.equilatero
        angLcdg = obj.P(3) + pi + [pi/6 -7*pi/6 -pi/2]; 
    else
        if obj.nLados == 3
            % Variables locales de la función
            mitheta = obj.theta;
            migamma = obj.gamma;
            angLcdg = obj.P(3) + [pi + migamma(1), ...
                       - mitheta(2) + migamma(2),...
                       pi - mitheta(2) - mitheta(3) + migamma(3)];
        end
    end
end

%% Parámetro B
function B = get.B(obj)
    % Coordenadas globales de los vértices de la plataforma. Se establecen
    % con P
    if obj.equilatero
            phi = obj.P(3);     % Por comodidad, llamamos phi a la orientación
            
            B(1,1) = obj.P(1) - obj.L(1,3)*(cos(pi/6 + phi));
            B(1,2) = obj.P(2) - obj.L(1,3)*(sin(pi/6 + phi));

            B(2,1) = B(1,1) + obj.T(2)*cos(phi);
            B(2,2) = B(1,2) + obj.T(2)*sin(phi);

            B(3,1) = B(1,1) + obj.T(3)*cos(pi/3 + phi);
            B(3,2) = B(1,2) + obj.T(3)*sin(pi/3 + phi);
    else      
        if obj.nLados == 3
            B = obj.P(1:2) + obj.L .* [cos(obj.angLcdg); sin(obj.angLcdg)]';
        end
    end
end


%% Asignación de objetos dependientes para indicar el error
% Las propiedades dependientes no se pueden asignar, pero las siguientes
% funciones permiten personalizar el mensaje de error cuando se trata de
% realizar esta asignación
function obj = set.nLados(~, ~)  
    error('nLados no se puede asignar');            
end
function obj = set.theta(~, ~)  
    error('theta no se puede asignar');            
end
function obj = set.coord(~, ~)  
    error('coord no se puede asignar');            
end
function obj = set.cdg(~, ~)  
    error('cdg no se puede asignar');            
end
function obj = set.L(~, ~)  
    error('L no se puede asignar');            
end
function obj = set.gamma(~, ~)
    error('gamma no se puede asignar');            
end
function obj = set.angLcdg(~, ~)
    error('angLcdg no se puede asignar');            
end
function obj = set.B(~, ~)
    error('B no se puede asignar');            
end


%% Métodos extras para manejar y dibujar el robot


%% Dibujo de la plataforma en ejes locales
function plotLocal(obj)    
    % Dibujo de la plataforma en ejes locales
    figure
    if obj.nLados == 3
        x = [obj.coord(:,1)]; 
        y = [obj.coord(:,2)];
        fill(x, y, [0.8 0.8 0.8]);
        hold on
        plot(obj.cdg(1), obj.cdg(2), 'rx');
        axis equal
    end
end


%% Dibujo de la plataforma en ejes globales
% function plot(obj)    
%     % Dibujo de la plataforma en ejes globales
%     figure
%     if obj.nLados == 3
%         x = [obj.B(:,1)];
%         y = [obj.B(:,2)];
%         fill(x, y, [0.8 0.8 0.8]);
%         hold on
%         plot(obj.P(1), obj.P(2), 'rx');
%         axis equal
%     end
% end


end             % Fin de los métodos
end             % Fin de la clase