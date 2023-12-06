function aSolve = glauertSolution(K, a)

% Define una tolerancia para la convergencia
tol = 1e-6;

% Define un valor inicial para "a" (puede ser cualquier valor inicial razonable mayor a 1/3)
a0  = a;

% Inicializa la diferencia entre las iteraciones sucesivas
dif = Inf;

% Define la ecuación que deseas resolver
F  = @(a) (3/4)*K * a^3 - (1 + (5/4) * K) * a^2 + a * (K + 2) - 1;
% p  = [(3/4)^K,  -(1 + (5/4) * K), (K + 2), - 1];
% r  = roots(p);

% Define la derivada de la ecuación con respecto a "a"
dF = @(a) 3 * (3/4)*K * a^2 - 2 * (1 + (5/4) * K) * a + K + 2;

% Aplica el método de Newton-Raphson
while dif > tol
    aSolve = a0 - F(a0) / dF(a0);
    dif    = abs(aSolve - a0);
    a0     = aSolve;
end

% err = r - a;

end