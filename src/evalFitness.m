% Evalua el fitness de un número NRobots de robots
function [fitness, malos] = evalFitness(NRobots, NPadres, malos, niter, fitness, tipoNorm, qd, robot)
    
    % EN la primera iteración externa, se evalúan todos los robots. A
    % partir de la siguiente, solo los hijos
    persistent vectIter;
    vectIter = 1:NRobots;
    if niter == 2
        vectIter = NPadres:NRobots;
    end
    
    for k = vectIter
        % Actualización de los parámetros, para que se correspondan a esas
        % coordenadas
        if length(robot(k).P) > 3
            robot(k).sgn = robot(k).P(4:6);
        else
            signos = randi([1 8]);
            switch signos
                case 1
                    robot(k).P(4:6) = [1 1 1];
                case 2
                    robot(k).P(4:6) = [1 1 -1];
                case 3
                    robot(k).P(4:6) = [1 -1 1];
                case 4
                    robot(k).P(4:6) = [1 -1 -1];
                case 5
                    robot(k).P(4:6) = [-1 1 1];
                case 6
                    robot(k).P(4:6) = [-1 1 -1];
                case 7
                    robot(k).P(4:6) = [-1 -1 1];
                case 8
                    robot(k).P(4:6) = [-1 -1 -1];
            end
            robot(k).sgn = robot(k).P(4:6);
        end
        
        % Los padres ya tienen el MCI hecho de la ronda anterior
        
        % Calculamos el MCI
        robot(k) = MCI(robot(k));
        
        % Si algún ángulo es complejo, borramos el robot y generamos otro
        if sum(isnan(robot(k).q)) > 0
            robot(k) = generaRobots(1, robot(k), robot(k).P(1:3));
            malos = malos + 1;
        end
        
         % Inverso del error cuadrático entre los q producto del MCI y los q dato
        fitness(k,2) = 1/norm(qd - robot(k).q, tipoNorm);

    end
end