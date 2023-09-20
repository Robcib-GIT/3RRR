% Genera Num robots válidos, con MCI incluído
function robotGen = generaRobots(Num, modelo, P)
    % Configuración
    L = modelo.L;
    T = modelo.T;
    O = modelo.O;
    area = abs(cross([O(1,:) - O(2,:), 0], [O(1,:) - O(3,:), 0]));
    area = area(3);
    malos = 0;

    k = Num;
    while k
        % Generamos Num robots 
        robotGen(k) = robot3RRR;
        robotGen(k).L = L;
        robotGen(k).T = T;
        robotGen(k).O = O;

        % Generación de x, y, theta
        if malos >= 2* Num
            robotGen(k).P(1) = randBtw(0.075, 0.44);           % Intervalo [0.09 0.5]
            robotGen(k).P(2) = randBtw(-0.015, 0.36);          % Intervalo [0 0.41]
            robotGen(k).P(3) = randBtw(-2*pi/3, 2*pi/3);       % Intervalo +/- 2*pi/3
        else
            robotGen(k).P = P + 0.1*area*randn(1,3);
        end
        
        % Generación de los signos
        signos = randi([1 8]);
        switch signos
            case 1
                robotGen(k).P(4:6) = [1 1 1];
            case 2
                robotGen(k).P(4:6) = [1 1 -1];
            case 3
                robotGen(k).P(4:6) = [1 -1 1];
            case 4
                robotGen(k).P(4:6) = [1 -1 -1];
            case 5
                robotGen(k).P(4:6) = [-1 1 1];
            case 6
                robotGen(k).P(4:6) = [-1 1 -1];
            case 7
                robotGen(k).P(4:6) = [-1 -1 1];
            case 8
                robotGen(k).P(4:6) = [-1 -1 -1];
        end  

        % Calculamos el MCI
        robotGen(k) = MCI(robotGen(k));

        % Si algún ángulo es complejo, borramos el robotGen y generamos otro
        if sum(isnan(robotGen(k).q)) > 0
            malos = malos + 1;
            continue;
        end

        % Decrementamos
        k = k - 1;

    end
end

