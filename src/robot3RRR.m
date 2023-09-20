%% Robot 3RRR 
% 
% TFG de Jorge García Samartín
%
%% Esquema del robot
% 
%  
% <<../figuras/esquema.png>>
% 
classdef robot3RRR < handle
% Clase para gestionar el robot
% -------------------------------------------------------------------------
% J.G.S. v1.01 2020-01-30
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Declaración de propiedades visibles accesibles desde fuera de la clase
% -------------------------------------------------------------------------
properties
    O = [0 0; 0.595 0; 0.297 0.515];% Coordenadas de los apoyos
    L = [0.2 0.2 0.1155;            % Longitud de la barra j-ésima de la articulación i.           
         0.2 0.2 0.1155;
         0.2 0.2 0.1155];           % Longitudes de los barras hasta P
    T = [0.2 0.2 0.2];              % Longitudes de los lados del triángulo central

    q = [1.529 1.588 -2.604];       % Angulos de la primera barra respecto horizontal
    P = [0.3 0.2 0];                % Coordenadas del punto central
    sgn = [1 1 1];
    psi = zeros(3,1);
    mierror = zeros(3,1);
    dib = 0;                        % Si queremos que se dibuje o no el robot
    graficaRobot                    % Axes de la figura
    drawW = 0;                      % 1 si queremos que pinte el workspace
    drawS = 0;                      % 1 si queremos que pinte los puntos singulares
    drawM = 0;                      % 1 si queremos que pinte un mapa de color, en función de cond(J)
    conv = 0;                       % 1 si converge el MCD
    fitMax = 0;                     % Mejor resultado (del genético)
    net = cell(8,1);                         % Redes neuronales
end

% -------------------------------------------------------------------------
% Declaración de propiedades dependientes que pueden asignarse
% -------------------------------------------------------------------------
properties (SetAccess = public, Dependent = true)
    A = zeros(3,2); 
    Plataforma
end
% -------------------------------------------------------------------------
% Declaración de propiedades ocultas accesibles desde fuera de la clase
% -------------------------------------------------------------------------
properties (Hidden, Constant)
    g  = 9.81;          % gravedad                                  (m/s^2)
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
    hayCambios
end


%% Declaración de métodos
methods

%% Constructor. 
% Método principal de la clase. Se ejecuta cuando se carga un nuevo objeto
%
% Puede llamarse mediante;
%   obj = robot3RRR(O,L,T,q,P,s,sgn)
%
function obj = robot3RRR(varargin) %                 (L, T, O, s, sgn, error, q, P, dib)
    
    % Si se han pasado datos al constructor
    if nargin
        % Hacer un parse a los datos de entrada
        % http://es.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
        p = inputParser;

        % Declarar parámetros opcionales en la creación del objeto
        addOptional(p, 'O',   [], @isnumeric);
        addOptional(p, 'L',   [], @isnumeric);
        addOptional(p, 'T',   [], @isnumeric);
        addOptional(p, 'q',   [], @isnumeric);
        addOptional(p, 'P',   [], @isnumeric);
        addOptional(p, 'sgn', [], @isnumeric);
        addOptional(p, 'dib', [], @isnumeric);
        addOptional(p, 'dW',  [], @isnumeric);
        addOptional(p, 'dS',  [], @isnumeric);
        addOptional(p, 'dM',  [], @isnumeric);

        % Procesa los parámetros de entrada
        p.KeepUnmatched = false;
        parse(p, varargin{:})

        % Asignar las propiedades del objeto
        obj.O =  p.Results.O;
        obj.L =  p.Results.L;
        obj.T =  p.Results.T;
        obj.q =  p.Results.q;
        obj.P =  p.Results.P;
        obj.dib =  p.Results.dib;
        obj.sgn =  p.Results.sgn;
        obj.drawW  =  p.Results.dW;
        obj.drawS  =  p.Results.dS;
        obj.drawM  =  p.Results.dM;

    end

    % Redes Neuronales
    load('./nets/net_11-1.mat');
    obj.net{1} = net;
    load('./nets/net_1-11.mat');
    obj.net{2} = net;
    load('./nets/net_1-1-1.mat');
    obj.net{3} = net;
    load('./nets/net_111.mat');
    obj.net{4} = net;
    load('./nets/net_-111.mat');
    obj.net{5} = net;
    load('./nets/net_-11-1.mat');
    obj.net{6} = net;
    load('./nets/net_-1-11.mat');
    obj.net{7} = net;
    load('./nets/net_-1-1-1.mat');
    obj.net{8} = net;
    
    % Configuración de la figura y primer gráfico
    if obj.dib
        f = figure;
        obj.graficaRobot = axes(f);

        axis(obj.graficaRobot, 'equal');
        obj.robotPlot;
    end
    
    % Añadimos un listener para redibujar el robot cuando haya cambios
    if obj.dib
        addlistener(obj,'hayCambios',@(~,~) obj.robotPlot);
    end

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

%% Parámetro O. Coordenadas de los apoyos
function obj = set.O(obj, O)
    if isempty(O) 
       obj.O = obj.O;
    else
       obj.O =  O;
    end
    notify(obj,'hayCambios');
end

function O = get.O(obj)
    O = obj.O;
end

%% Parámetro L. Longitud de los brazos desde apoyos al centro P
function obj = set.L(obj, L)
    if isempty(L)
        obj.L = obj.L;
    else
        obj.L = L;
    end
    
    obj.L(:,3) = obj.Plataforma.L;
    
    notify(obj,'hayCambios');
end

function L = get.L(obj)
    L = obj.L;
end

%% Parámetro T. 
function obj = set.T(obj, T)
    % Longitud de los lados del triángulo central
    if isempty(T)
        obj.T = obj.T;
    else
        % Comprobar que la suma de 2 lados es mayor que el tercero
%         if sum(T(1:2) > T(3)) || sum(T(2:3) > T(1)) || sum(T(1:2:3) > T(2))
%             error('Tríangulo central mal definido')
%         else      
%             obj.T = T;
%         end
        obj.T = T;
    end
    notify(obj,'hayCambios');
end

function T = get.T(obj)
    T = obj.T;
end

%% Parámetro q 
function obj = set.q(obj, q)
    if isempty(q)
        obj.q = [1.529 1.588 -2.604];
        %obj.q =  input('q: ');
    else
        obj.q = q;
    end
    notify(obj,'hayCambios');
end

function q = get.q(obj)
    q = obj.q;
end

%% Parámetro P 
function obj = set.P(obj, P)
    if isempty(P)
        obj.P = obj.P;
    else
        obj.P = P;
    end
    notify(obj,'hayCambios');
end

function P = get.P(obj)
    P = obj.P;
end

%% Parámetro s
% function obj = set.s(obj, s)
%     if isempty(s)
%         obj.s = obj.s;
%     else
%         obj.s = s;
%     end
%     notify(obj,'hayCambios');
% end
% 
% function s = get.s(obj)
%     s = obj.s;
% end

%% Parámetro sgn
function obj = set.sgn(obj, sgn)
    if isempty(sgn)
        obj.sgn = obj.sgn;
    else
        obj.sgn = sgn;
    end
    notify(obj,'hayCambios');
end

function sgn = get.sgn(obj)
    sgn = obj.sgn;
end

%% Parámetro dib
function obj = set.dib(obj, dib)
    if isempty(dib)
        obj.dib = obj.dib;
    else
        obj.dib = dib;
    end
    notify(obj,'hayCambios');
end

function dib = get.dib(obj)
    dib = obj.dib;
end

%% Ángulos psi (solo se ponen a cero para que no salte error)
function obj = set.psi(obj, psi)  
    if isempty(psi)
        obj.psi = obj.psi;
    else
        obj.psi = psi;
    end
    %notify(obj,'hayCambios');           
end

function psi = get.psi(obj)
    psi = obj.psi;
end

%% Error
function obj = set.mierror(obj, mierror)
    if isempty(mierror)
        obj.mierror = obj.mierror;
    else
        obj.mierror = mierror;
    end
end

function mierror = get.mierror(obj)
    mierror = obj.mierror;
end
%% Ejes de la figura
function obj = set.graficaRobot(obj, graficaRobot)
    if isempty(graficaRobot)
        obj.graficaRobot = obj.graficaRobot;
    else
        obj.graficaRobot = graficaRobot;
    end
end

function graficaRobot = get.graficaRobot(obj)
    graficaRobot = obj.graficaRobot;
end
%% Dibujar workspace
function obj = set.drawW(obj, drawW)
    if isempty(drawW)
        obj.drawW = obj.drawW;
    else
        obj.drawW = drawW;
    end
end

function drawW = get.drawW(obj)
    drawW = obj.drawW;
end
%% Dibujar singularidades
function obj = set.drawS(obj, drawS)
    if isempty(drawS)
        obj.drawS = obj.drawS;
    else
        obj.drawS = drawS;
    end
end

function drawS = get.drawS(obj)
    drawS = obj.drawS;
end

%% Propiedades dependientes (solo tienen método get)
% Las propiedades que se definen a continuación son dependientes de los
% datos de entrada y se calculan cada vez que cambiamos uno de ellos. Por
% tanto no se les puede asignar un valor como a los datos de entrada. Ello
% implica que son de solo de lectura, tienen método GET pero no método SET

%% Parámetro A. Posiciones de las articulaciones
function A = get.A(obj)
     miq = [obj.q; obj.q - pi/2];
     A = obj.O + obj.L(:,1) .* cos(miq');
end

%% Parámetro Platform- Objeto con la plataforma
function Plataforma = get.Plataforma(obj)
    Plataforma = Platform(obj.T, obj.P);
end

%% Asignación de objetos dependientes para indicar el error
% Las propiedades dependientes no se pueden asignar, pero las siguientes
% funciones permiten personalizar el mensaje de error cuando se trata de
% realizar esta asignación
function obj = set.A(~, ~)  
    error('A no se puede asignar');            
end
function obj = set.Plataforma(~, ~)  
    error('Plataforma no se puede asignar');            
end

%% Métodos extras para manejar y dibujar el robot


%% Modelo cinemático inverso
% Cardona, Dimensional Synthesis
function obj = MCI(obj)
    
    % Variables locales
    miO = obj.O;
    miL = obj.L;
    miB = obj.Plataforma.B;
    misigno = obj.sgn;
    
    mipsi = zeros(3,1);
    
    % Inicialización de variables
    fracc = zeros(3, 1);    % Vector con los cosenos de los ángulos pasivos
    
    % Cerramos las ventanas abiertas
    %close all;

    for i = 1:3
        
        ix = miB(i,1) - miO(i,1);
        iy = miB(i,2) - miO(i,2);

        % Teorema del coseno
        fracc(i) = (ix^2 + iy^2 - miL(i,1)^2 - miL(i,2)^2) / (2*miL(i,1)*miL(i,2));
        mipsi(i) = misigno(i)*acos(fracc(i));
         
        % Psi tiene distinto criterio de signos en el tercer ángulo
%         if i == 3
%             mipsi(i) = -mipsi(i);
%         end
         
        % Parámetros de la ecuación cuadrática
        a = miL(i,1)^2 + miL(i,2)^2 + 2*miL(i,1)*miL(i,2)*fracc(i);
        b = -2*ix*(miL(i,1) + miL(i,2)*fracc(i));
        c = ix^2 - miL(i,2)^2*(1-fracc(i)^2);
        
        % Discriminamos la solución en función del signo de psi
        if misigno(i) > 0
            sol = acos((-b + sqrt(b^2 - 4*a*c))/2/a);
        else
            sol = acos((-b - sqrt(b^2 - 4*a*c))/2/a);
        end
        
        % Da error si la solución no es real
        if imag(sol) ~= 0
%             disp('Error');
            obj.psi(i) = NaN;
            obj.q(i) = NaN;
            obj.mierror(i) = 1;
            break;
        else     
            % Vemos si la solución tiene el signo del coseno correcto
            miA = miO(i,:) + [miL(i,1)*cos(sol) miL(i,1)*sin(sol)];
            
            if (abs(norm(miA - miB(i,:)) - miL(i,2)) < 1e-5)
                obj.q(i) = sol;
            else
                obj.q(i) = -sol;
            end
            obj.psi = mipsi;
        end
        
    end
end

%% Modelo cinemático directo
function obj = MCD(obj, varargin)  

    %% Parámetros de configuración
    % Si se han pasado datos al método
    if nargin > 1
        % Hacer un parse a los datos de entrada
        % http://es.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
        p = inputParser;

        % Declarar parámetros opcionales en la creación del objeto
        addOptional(p, 'MAX_ITER',  [], @isnumeric);
        addOptional(p, 'TOL',       [], @isnumeric);
        addOptional(p, 'TOL_SAL',   [], @isnumeric);
        addOptional(p, 'DRAWING',   [], @islogical);

        % Procesa los parámetros de entrada
        p.KeepUnmatched = false;
        parse(p, varargin{:})
        
        % Asignar las propiedades del objeto
        MAX_ITER    =  p.Results.MAX_ITER;  
        TOL         =  p.Results.TOL;  
        TOL_SAL     =  p.Results.TOL_SAL; 
        DRAWING     =  p.Results.DRAWING;

    else
        % Valores por defecto
        MAX_ITER = 100;     % Número máximo de iteraciones
        TOL = 1e-4;         % Error máximo permitido
        TOL_SAL = 1e-3;     % Error máximo que provoca error tras las MAX_ITER iteraciones
        DRAWING = 0;        % 1 si queremos que saque gráficas de convergencia
    end

    % Inicialización de variables
    F = zeros(3,1);     % Vector con las ecuaciones de los bucles cerrados
    Z = ones(3, 1);     % Vector con las incógnitas iteración nueva
    
    % Vector con los ángulos alpha
    %a = [pi/6 5*pi/6 3*pi/2];
    a = wrapToPi(obj.Plataforma.angLcdg + pi);
    
    % Método de Gauss-Newton
    
    % Valor inicial de Z
    Zn = obj.P';
    
    if DRAWING 
        u = zeros(MAX_ITER, 3);
    end
    
    % Cerramos todo
    close all;
    
    % Para un máximo de iteraciones (o hasta que el error sea suf
    % pequeño)
    niter = 1;
    while norm(Z-Zn) > TOL && niter < MAX_ITER
        
        % Actualizamos Z con el valor de la iteración anterior
        Z = Zn;
        
        fprintf('********************* ITERACION %2d ************************ \n', niter);
        % Vector de cadenas cerradas
        F = evalF(obj, Z, a);

        % Jacobiana
        J = evalJ(obj, Z, a);
        
        % Actualizamos el valor de las incógnitas
        Zn = Z - J \ F;
        
        % Para que phi no se desmadre
        Zn(3) = wrapToPi(Zn(3));
        
        if DRAWING
            u(niter,:) = Z';
        end
        
        niter =  niter + 1;        
    end
    
    % Cambiamos la posición del punto P 
    if imag(Zn) == 0
        obj.P = Zn';
    end
    
    if DRAWING
        u = u(1:niter-1, :); 
        
        % Si solo hay una interación, se hace una "chapuza" paa que salga
        % bien el gráfico
        if size(u, 1) == 1
            u(2,:) = u(1,:);
        end
        figure;
        plot(u);
        
        % Que solo haya enteros en el eje horizontal
        xt = xticks;
        xticks(floor(xt));
        
        % Formato
        title('Evolución de x, y, \theta');
        legend('x', 'y', '\theta');
    end
    
    % Comprobación de salida por exceso número de iteraciones
    if  niter >= MAX_ITER && norm(Z-Zn) > TOL_SAL
        %error('No ha convergido el modelo directo');
        obj.conv = obj.conv + 1;
        disp('No ha convergido el algortimo de MCD');
    else
        notify(obj,'hayCambios');
        obj.conv = 0;
    end
    
end

%% Funciones adicionales para el modelo cinemático directo
% Vector de cadenas cerradas
function F = evalF(obj, Z, a)

    % Por comodidad, damos nombre a las componentes de Z
    x = Z(1);
    y = Z(2);
    phi = Z(3);
    
    % Declaramos F por velocidad
    F = zeros(3, 1);
    
    for i = 1:3
        F(i) = (x - obj.O(i,1) - obj.L(i,1)*cos(obj.q(i)) - obj.L(i,3)*cos(a(i) + phi))^2 + ...
               (y - obj.O(i,2) - obj.L(i,1)*sin(obj.q(i)) - obj.L(i,3)*sin(a(i) + phi))^2 - obj.L(i,2)^2;
    end
    
end

% Jacobiana
function J = evalJ(obj, Z, a)

    % Por comodidad, damos nombre a las componentes de Z
    x = Z(1);
    y = Z(2);
    phi = Z(3);
    
    % Declaramos J para ganar velocidad
    J = zeros(3,3);
    
    for i = 1:3
        J(i,1) = 2*(x - obj.O(i,1) - obj.L(i,1)*cos(obj.q(i)) - obj.L(i,3)*cos(a(i) + phi));
        J(i,2) = 2*(y - obj.O(i,2) - obj.L(i,1)*sin(obj.q(i)) - obj.L(i,3)*sin(a(i) + phi));
        J(i,3) = J(i,1)*obj.L(i,3)*sin(a(i) + phi) - obj.L(1,3)*J(i,2)*cos(a(i) + phi);
    end
end

%% Modelo cinemático directo por algoritmos genéticos
% Parámetros
% NRobots: número de robots en cada generación
% PPadres: tanto por uno de robots escogidos como padres
% tipoNorm: norma vectorial que evalúa la función objetivo
% PMutar: probabilidad (0, 1) de que un robot mute
% PMigrar: probabilidad (0,1) de que un robot sea sustituído por otro
% aleatorio
% max_niter: máximo número de generaciones que se simulan
% tol: error máximo permitido para finalizar prematuramente la simulación
% tolSal: valor máximo de error por encima del cual se considera que el
% algoritmo no ha convergido, aún pasadas las max_niter iteraciones
% dibujar: vector de cinco elementos que indican, cuando están a 1, que se
% desea dibujar histograma, evolución de las soluciones, perceniles, % robots
% defectuosos y la distribución final de los mismos
function obj = MCDGen(obj, varargin)

    %% Lectura de los argumentos
    % Si se han pasado datos al método
    if nargin > 1
        % Hacer un parse a los datos de entrada
        % http://es.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
        p = inputParser;

        % Declarar parámetros opcionales en la creación del objeto
        addOptional(p, 'NRobots',   [], @isnumeric);
        addOptional(p, 'PPadres',   [], @isnumeric);
        addOptional(p, 'tipoNorm',  [], @isnumeric);
        addOptional(p, 'PMutar',    [], @isnumeric);
        addOptional(p, 'PMigrar',   [], @isnumeric);
        addOptional(p, 'max_niter', [], @isnumeric);
        addOptional(p, 'tol',       [], @isnumeric);
        addOptional(p, 'tolSal',    [], @isnumeric);
        addOptional(p, 'dibujar',   [], @isvector);

        % Procesa los parámetros de entrada
        p.KeepUnmatched = false;
        parse(p, varargin{:})
        
        % Asignar las propiedades del objeto
        NRobots   =  p.Results.NRobots;  
        PPadres   =  p.Results.PPadres; 
        tipoNorm  =  p.Results.tipoNorm;
        PMutar    =  p.Results.PMutar;  
        PMigrar   =  p.Results.PMigrar;  
        max_niter =  p.Results.max_niter;
        tol       =  p.Results.tol;      
        tolSal    =  p.Results.tolSal;   
        dibujar   =  p.Results.dibujar; 

    else
        % Valores por defecto
        NRobots   = 75;
        PPadres   = 0.2;
        tipoNorm  = inf;
        PMutar    = 0.7;
        PMigrar   = 0;
        max_niter = 100;
        tol       = 5e-3;           % Resolución del encoder
        tolSal    = 17e-3;          % 1 grado
        dibujar   = [0 0 0 0 0];    % Dibujar solo la evolución
    end
    
    %% Parámetros de configuración
    
    % Valor deseado
    qd = obj.q;
    
    % Comprobación inicial de que no se está llevando al robot al punto en
    % el que se encuentra
    obj.MCI;
    if obj.q == qd
        return;
    end
    
    % Restauramos qd
    obj.q = qd;
    
    % Copias de seguridad del dibujo
    dibCopia = obj.dib;
    grafCopia = obj.graficaRobot;

    % Vector con los errores cuadráticos de cada robot
    fitness = zeros(NRobots, 2);  
    % La primera columna de fitness es ell id del robot (para poder ordenar)
    fitness(:,1) = 1:NRobots;

    % Para el dibujado
    controlAc = zeros(max_niter, 4);
    percentil = zeros(max_niter, 4);

    % Cada cuantas interaciones dibujar la distribución de robots
    ciclosDistr = 0;

    % Número total de robots malos
    malosTot = 0;
    
    % Conversión de porcentaje a número de padres
    NPadres = floor(NRobots * PPadres);
    
    % Definición de los parámetros de dibujado
    dibHistograma = dibujar(1);
    dibEvolucion = dibujar(2);
    dibPercentiles = dibujar(3);
    dibMalos = dibujar(4);
    dibDistrFinal = dibujar(5);
    printResults = 0;
    
    % La primera columna de fitness contiene los id de los robots
    fitness(:,1) = 1:NRobots;
    
    % Cerramos todo
    close all;
    
    conv = 1;

    %% Generación
    robot = generaRobots(NRobots, obj, obj.P(1:3));

    for niter = 1:max_niter
        fprintf('********************* ITERACION %2d ************************ \n', niter);

        % Antes de evaluar, ponemos malos a cero
        malos = 0;

        %% Evaluación
        [fitness, malos] = evalFitness(NRobots, NPadres, malos, niter, fitness, tipoNorm, qd, robot); 

        % Mejor robot
        M = max(fitness(:,2));

        % Sumamos los malos
        malosTot = malosTot + malos;

        % Parámetros del dibujo
        if dibEvolucion
            m = min(fitness(:,2));
            med = mean(fitness(:,2));
            controlAc(niter,:) = [M m med malos];
        end

        if (M >= 1/tol)
            break;
        end

        %% Selección

        % Elegimos a los NPadres mejores como padres
        fitOrd = sortrows(fitness, 2, 'descend');

        %clear padre;
        padre = robot(fitOrd(1:NPadres,1));

        % Esos padres pasan a la siguiente generación
        clear robot;
        robot(1:NPadres) = padre;

        % Para dibujar los percentiles
        if dibPercentiles
            for i = 1:4
                mejoresP = maxk(fitness(:,2), floor(NPadres*i/4));
                percentil(niter,i) = fitness(find(fitness(:,2) == mejoresP(i), 1),2);
            end
        end
        
        % Para no evaluar a los padres de nuevo
        fitness(1:NPadres,2) = fitOrd(1:NPadres,2);

        %% Cruce, mutación y migración

        % Cruce por media aritmética (P) y en un punto (sgn)
        for i = NPadres:NRobots

            % Seleccionamos los padres
            p1 = randi([1 NPadres]);

            % Es una especie de do-while
            p2 = p1;
            while p2 == p1
                p2 = randi([1 NPadres]);
            end

            ptCruce = randi([4 5]); 
            robot(i).P(1:3) = 0.5 * (padre(p1).P(1:3) + padre(p2).P(1:3));
            robot(i).P(4:6) = [padre(p1).P(4:ptCruce) padre(p2).P(ptCruce + 1:end)];
            
            % Mutación
            if rand < PMutar
                robot(i).P(1:2) = robot(i).P(1:2) .* randBtw(.9, 1.1, 1, 2);
                robot(i).P(3) = robot(i).P(3) .* randBtw(.9, 1.1);
            end
            
            % Migración
            if rand < PMigrar
                robot(i) = generaRobots(1, obj, obj.P(1:3));
            end            

        end
        %% Dibujado de la distribución de robots
        if ciclosDistr && ~mod(niter, ciclosDistr)
            figure;
            hold on;
            for i=1:NRobots
                plot(robot(i).P(1), robot(i).P(2), 'x');
            end
        end    

    end    
        
    % Mejor valor encontrado de fitness
    fitMax = M^-1;

    % Índice del robot de mejor fitness
    best = find(fitness(:,2) == M, 1);

    % En P solo guardamos las componentes reales (no los signos)
    clear obj.P;
    obj.P = robot(best).P(1:3);
    
    %obj
    
%     obj.conv = conv;
%     obj.it = niter;
    
    % Restauramos las copias de seguridad
    obj.dib = dibCopia;
    obj.graficaRobot = grafCopia;
    
    %obj
    
    % Comprobación de salida por exceso número de iteraciones
    if  niter >= max_niter && fitMax > tolSal
        robot(best).P
        obj.fitMax = fitMax
        obj.conv = 0;
        %error('No ha convergido el modelo directo por algoritmos genéticos');
    else
        obj.fitMax = fitMax
        obj.conv = 1;
    end
    
    % Notificamos de que hay cambios (para el dibujado)
    notify(obj,'hayCambios');
        
    if printResults
        mejQ = robot(best).q
        mejP = robot(best).P
    end

    if dibHistograma
        histogram(fitness(:,2), 10);
        title('Histograma');
        ylabel('Número de robots');
        xlabel('Inverso del fitness');
    end

    if dibEvolucion
        figure;
        hold on;
        plot(controlAc(1:niter,1));
        plot(controlAc(1:niter,2));
        plot(controlAc(1:niter,3));
        xlim([1 niter-1]);
        title('Evolución');
        legend('Máximo', 'Mínimo', 'Media', 'Location', 'NorthWest');
        ylabel('Inverso del fitness');
        xlabel('Número de generaciones');
    end

    if dibMalos
        figure;
        hold on;
        plot(100*controlAc(2:niter,4)/NRobots);
        title('Robots defectuosos');
        legend('%Malos', 'Location', 'southeast');
        ylabel('% Robots defectuosos');
        xlabel('Número de generaciones');
    end

    if dibPercentiles
        figure;
        hold on;
        plot(percentil(1:niter,1));
        plot(percentil(1:niter,2));
        plot(percentil(1:niter,3));
        plot(percentil(1:niter,4));
        xlim([1 niter-1]);
        title('Evolución del fitness de los padres');
        legend('Cuartil 1', 'Cuartil 2', 'Cuartil 3', 'Cuartil 4', 'Location', 'NorthWest');
        ylabel('Inverso del fitness');
        xlabel('Número de generaciones');
    end

    if dibDistrFinal
        figure;
        hold on;
        for i=1:NRobots
            plot(robot(i).P(1), robot(i).P(2), 'dr');
        end
        title('Posición (x e y) de los robots de la última generación');
        ylabel('Coordenada y');
        xlabel('Coordenada x');
        
    end

end

%% Modelo cinemático directo por redes neuronales
function obj = MCDNeur(obj, varargin)
    
    Puntos = zeros(8,3);
    for i = 1:8
        Puntos(i,1:3) = sim(obj.net{i}, obj.q');
    end

    % Medimos distancia en P
    for j = 1:length(Puntos)
        Puntos(j,4) = norm(obj.P(1:3) - Puntos(j,1:3));
    end
    min1 = find(Puntos(:,4) == nanmin(Puntos(:,4)));
    obj.P = Puntos(min1(1),1:3);

   % Signo
   conf = [1 1 -1;
        1 -1 1;
        1 -1 -1;
        1 1 1;
        -1 1 1;
        -1 1 -1;
        -1 -1 1;
        -1 -1 -1];
   obj.sgn = conf(min1,:);
   
end
%% Mide en puntos (con resolución res) el tamaño del espacio de trabajo del
% robot
function [workspace, params, desvs, sings, puntos, lista] = medirParams(obj, phi, sgn, varargin)
    
    %close all;
    
    % Copia de seguridad
    copiaP = obj.P;
    copiaQ = obj.q;
    copiasgn = obj.sgn;
    
    % Configuración del robot
    obj.P(3) = phi;
    obj.sgn = sgn;   
    
    % Parámetros de configuración de la medida
        % Hacer un parse a los datos de entrada
    % http://es.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
    p = inputParser;

    % Declarar parámetros opcionales en la creación del objeto
    addOptional(p, 'res',        [], @isnumeric);
    addOptional(p, 'margen',     [], @isnumeric);
    addOptional(p, 'porcManip',  [], @isnumeric);
    addOptional(p, 'tolP',       [], @isnumeric);
    addOptional(p, 'tolQ',       [], @isnumeric);
    addOptional(p, 'ejesGraf',   []);

    % Procesa los parámetros de entrada
    p.KeepUnmatched = false;
    parse(p, varargin{:})

    % Asignar las propiedades del objeto
    % Resolución de la medida (.08 para que se vea lleno)
    if p.Results.res
        res       =  p.Results.res;
    else
        res = 0.015;
    end
    % Margen de la zona de medida respecto a las bases
    if p.Results.margen
        margen    =  p.Results.margen;
    else
        margen = 0.15;
    end
    % Tanto por uno de puntos DEL WORKSPACE que queremos medir
    if p.Results.porcManip
        porcManip =  p.Results.porcManip;
    else
        porcManip = 0.9;
    end
    % Valor por debajo del cual el determinate de Jp se considera cero
    if p.Results.tolP
        tolP      =  p.Results.tolP; 
    else
        tolP = 1e-4;
    end
    % Valor por debajo del cual el determinate de Jq se considera cero
    if p.Results.tolQ
        tolQ      =  p.Results.tolQ;   
    else
        tolQ = 1e-5;
    end
    % Ejes en los que se dibuja
    if isgraphics(p.Results.ejesGraf)
        ejesGraf     =  p.Results.ejesGraf; 
        hayEjes = 1;
    else
        ejesGraf = [];
        hayEjes = 0;
    end
    
    % Límites de la medida
    ejeX = [obj.O(1,1) obj.O(2,1)] + [-margen margen]; 
    ejeY = [obj.O(1,2) obj.O(3,2)] + [-margen margen]; 
    
    % Vectores para los bucles
    bucleX = ejeX(1):res:ejeX(2);
    bucleY = ejeY(1):res:ejeY(2);
    
    % Número de puntos a medir
    area = length(bucleX) * length(bucleY);
    
    % Variables de cuenta   
    workspace = 0;          % Tamaño del workspace
    med = 0;                % Número de puntos del workspace medidos
    singP = 0;              % Número de singularidades directas
    singQ = 0;              % Número de singularidades inversas
    singJ = 0;              % Número de singularidades combinadas
    eta = zeros(area,1);    % Inverso del número de condición (destreza)
    mu = zeros(area,1);     % Manipulabilidad
    gamma = zeros(area,1);  % Ganancia
    puntos = [];            % Vector con todas las P y q del workspace
    
    % Gráfica
    if (obj.drawW || obj.drawS || obj.drawM) && ~hayEjes
        figure;
        ejesGraf = axes;
        hold on;
    end
    
    % Medidas
    for x = bucleX
        for y = bucleY
            
            % Se lleva al robot a esa posición
            obj.P(1:2) = [x y];
            obj.MCI;
            
            % Si está en el workspace
            if ~isnan(sum(obj.q))
                % Medimos el workspace
                workspace = workspace + 1;
                
                % Añadimos P y q al vector de puntos
                puntos = [puntos; obj.P obj.q 0 0 0];
                
                if obj.drawW
                    plot(ejesGraf, x, y, '.', 'Color', 200/255*[1 1 1]);
                end
                
                % Analizamos si hay puntos singulares
                [Jp, Jq] = evalJpq(obj);
                J = Jq\Jp;
                
                if abs(det(Jp)) < tolP
                    singP = singP + 1;
                    puntos(end,7) = abs(det(Jp));
                    if obj.drawS
                        plot(ejesGraf, x, y, '.r');
                    end
                end
                
                if abs(det(Jq)) < tolQ
                    singQ = singQ + 1;
                    puntos(end,8) = abs(det(Jq));
                    if obj.drawS
                        plot(ejesGraf, x, y, '.k');
                    end
                end
                
                if abs(det(J)) < tolP*tolQ
                    singJ = singJ + 1;
                    puntos(end,9) = abs(det(J));
                    if obj.drawS
                        plot(ejesGraf, x, y, '.g');
                    end
                end
                
                % Si toca, pintamos el mpaa de color
                if obj.drawM
                    plot(ejesGraf, x, y, '.', 'Color', [0 255 0]/255 + 1/cond(J)/255*[255 -255 255]);
                end                    
                
                % Medimos la precisión
                if rand() < porcManip
                    
                    % Contamos que se ha medido un punto más
                    med = med + 1;
                                 
                    % Parámetros
                    eta(med) =  cond(J)^-1;
                    mu(med) = det(J);
                    gamma(med) = norm([norm(J(:,1),1) norm(J(:,2),1) norm(J(:,3),1)]);
                end
            end
            
        end
    end
    
    % Asignación de los resultados
    % Destreza, manipulabilidad y ganancia
    params(1) = mean(eta(1:med));
    params(2) = mean(mu(1:med));
    params(3) = mean(gamma(1:med));
    
    % Desviaciones típicas
    desvs(1) = std(eta(1:med));
    desvs(2) = std(mu(1:med));
    desvs(3) = std(gamma(1:med));
    
    % Singularidades
    sings(1) = singJ;
    sings(2) = singP;
    sings(3) = singQ;
    
    % Si hay que dibujar
    if obj.drawW || obj.drawS || obj.drawM
        xlim(ejeX);
        ylim(ejeY);
        axis equal;
        grid;
    end
    
    % Dibujo de la base
    if obj.drawW || obj.drawS || obj.drawM
        for i = 1:3
            plot(ejesGraf, obj.O(i,1), obj.O(i,2), '^b', 'MarkerSize', 10);
        end
    end
    
    % Restauramos la copia
    obj.P = copiaP;
    obj.q = copiaQ;
    obj.sgn = copiasgn;

end

%% Dibujo del robot
function robotPlot(obj)    

    %% Borramos los ejes
    if isempty(obj.graficaRobot)
        f = figure;
        obj.graficaRobot = axes(f);

        axis(obj.graficaRobot, 'equal');
    end
%     obj.graficaRobot.cla;
%     set(obj.graficaRobot,'nextplot','replacechildren');
    
    %% Parámetros
    % Margen de la zona de medida respecto a las bases
    margen = 0.15;
    
    % Límites de los ejes
    ejeX = [obj.O(1,1) obj.O(2,1)] + [-margen margen]; 
    ejeY = [obj.O(1,2) obj.O(3,2)] + [-margen margen];  
    
    %% Workspace y singularidades
    % Si procede, se dibujan el workspace y/o las singularidades
    if obj.drawW || obj.drawS
        obj.medirParams(obj.P(3), obj.sgn, 'ejesGraf', obj.graficaRobot);
    end
    
    %% Dibujado
    
    % Extremo terminal
    plot(obj.graficaRobot, obj.P(1), obj.P(2), 'ro');
    
    % Hay que ponerlo aquí, si no, no funciona
    hold(obj.graficaRobot, 'on');

    % Triángulo central
    for i = 1:3
        plot(obj.graficaRobot, obj.Plataforma.B(i,1), obj.Plataforma.B(i,2), 'bo');
    end
    % Actuadores
    for i = 1:3
        plot(obj.graficaRobot, obj.O(i,1), obj.O(i,2), 'kx');
    end

    % Articulaciones
    for i = 1:3
        plot(obj.graficaRobot, obj.A(i,1), obj.A(i,2), 'bo');
    end

    % Barras
    for i = 1:3
        plot(obj.graficaRobot, [obj.O(i,1) obj.A(i,1)], [obj.O(i,2) obj.A(i,2)], 'b');
        plot(obj.graficaRobot, [obj.Plataforma.B(i,1) obj.A(i,1)], [obj.Plataforma.B(i,2) obj.A(i,2)], 'b');
    end
    
    plot(obj.graficaRobot, [obj.Plataforma.B(1,1) obj.Plataforma.B(2,1)], [obj.Plataforma.B(1,2) obj.Plataforma.B(2,2)], 'b');
    plot(obj.graficaRobot, [obj.Plataforma.B(2,1) obj.Plataforma.B(3,1)], [obj.Plataforma.B(2,2) obj.Plataforma.B(3,2)], 'b');
    plot(obj.graficaRobot, [obj.Plataforma.B(3,1) obj.Plataforma.B(1,1)], [obj.Plataforma.B(3,2) obj.Plataforma.B(1,2)], 'b');
    
    grid(obj.graficaRobot, 'on');
    hold(obj.graficaRobot, 'off');
    
    % Límites
    xlim(obj.graficaRobot, ejeX);
    ylim(obj.graficaRobot, ejeY);
    
end

end             % Fin de los métodos
end             % Fin de la clase