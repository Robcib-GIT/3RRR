% Devuelve un número aleatorio en un intervalo (a,b).
%
%   r = randBtw(a, b, m, n)
%   Parámetros:
%       (a,b): extremos del intervalo (puede darse un vector)
%       (m,n): dimensiones de la matriz que se devuelve. Si solo se
%       especifica m, se devuelve una matriz cuadrada mxm.
function r = randBtw(a, b, m, n)

    % Si el primer argumento es un vector, se considera que se da a y b
    if length(a) > 1
        b = a(end);
        a = a(1);
        extraParams = 1;
    else
        extraParams = 0;
    end
    
    % Como en rand, si no se especifican columnas, se supone matriz
    % cuadrada y si no se especifica nada, entero
    if nargin == (3 - extraParams)
        n = m;
    elseif nargin == (2 - extraParams)
        m = 1;
        n = 1;
    end
    
    r = a + (b-a)*rand(m, n);
end