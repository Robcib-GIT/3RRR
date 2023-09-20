res = 0.01;

r = robot3RRR;
r.L = [.16 .18 .1155;.16 .18 .1155;.16 .18 .1155];
r.O = [0 0; .5 0; .25 .43];
r.sgn = [1 -1 1];

ejeX = [r.O(1,1) r.O(2,1)] + [-0.2 0.2]; 
ejeY = [r.O(1,2) r.O(3,2)] + [-0.2 0.2]; 
    
bucleX = ejeX(1):res:ejeX(2);
bucleY = ejeY(1):res:ejeY(2);

toPlot = [0 0 0 0 0 0];

max_error = 1.5;
metodo = 2;

figure;
hold on;
% Medidas
s = 0;
for x = bucleX
    for y = bucleY

        % Se lleva al robot a esa posición
        Por = [x y];
        r.P(1:2) = Por;
        r.P(3) = 0;
        r.MCI;
        qor = r.q';

        % Si está en el workspace
        if ~isnan(sum(r.q))
            
            switch metodo
                case 1
                    s = s+1;
                    Pnet = sim(net, r.q');
                case 2
                    r.P = [0 0 0];
                    s = s+1;
                    r.MCD;
                    Pnet = r.P';
                case 3
                    r.P = [0 0 0];
                    r.MCDGen;
                    Pnet = r.P';
            end
                      
            errorAbs = norm(Pnet(1:2) - Por');
            error = min(errorAbs, max_error);
            vecColor = [0 1 0] + error/max_error * [1 -1 0];
            toPlot = [toPlot; x y vecColor error];
            plot(x, y, '.', 'Color', vecColor);
        end
    end
    progreso = x/length(bucleX);
    disp(progreso);
end

for i = 1:3
    plot(r.O(i,1), r.O(i,2), 'xb', 'MarkerSize', 10);
end

xlim(ejeX);
ylim(ejeY);
axis equal;
grid;
disp(mean(toPlot(:,end)));