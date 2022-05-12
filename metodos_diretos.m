load "685_bus.mat"

A = Problem.A;
n = rows(A);

[L, U, P] = lu(A);
hf = figure();
spy(A);
print(hf, "658_bus.png", "-dpng");
hf = figure();
spy(L);
print(hf, "685_bus_L.png", "-dpng");
hf = figure();
spy(U);
print(hf, "685_bus_U.png", "-dpng");

fill = 100 - (nnz(A) / (nnz(L) + nnz(U))) * 100;
printf("Taxa de preenchimento = %f\n", fill);

solution = ones(n, 1);
b = A * solution;
x = A \ b;
printf("Solução aproximada para sistema linear = %f\n", norm(x, 2));
printf("Solução para sistema linear = %f\n", norm(solution, 2));

relDist = norm(solution - x, inf) / norm(solution, inf);
printf("Distância relativa = %f\n", relDist);

luRelDist = norm(A - P * L * U, inf) / norm(A, inf);
printf("Distância relativa da decomposição = %f\n", luRelDist);

deltaB = A * (solution - x);
independentRelDist = norm(deltaB, inf) / norm(b, inf);
printf("Distância relativa dos termos independentes = %f\n", independentRelDist);

residualNorm = norm(b - A * x, inf);
printf("Norma do resíduo = %f\n", residualNorm);

K = cond(A);
printf("Número de condicionamento = %f\n", K);

