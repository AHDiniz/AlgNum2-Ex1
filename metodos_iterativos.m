load "rail_1357.mat"

A = Problem.A;
n = rows(A);

solution = ones(n, 1);
b = A * solution;

isDiagDom = diagonal_dominante(A);
printf("É diagonal dominante? %s\n", isDiagDom);

tol = 0.0001;
printf("Tolerância = %f\n", tol);
nMaxIter = 1000;
printf("Número máximo de iterações = %f\n", nMaxIter);
omega = 1.5;
printf("Omega = %f\n", omega);

[xJacobi, erJacobi, iterJacobi, timeJacobi] = fast_jacobi(A, b, tol, nMaxIter);
normJacobi = norm(xJacobi, inf);
printf("Norma da solução (Jacobi) = %f\n", normJacobi);
printf("Erro (Jacobi) = %f\n", erJacobi(iterJacobi));
printf("Número de iterações (Jacobi) = %d\n", iterJacobi);
printf("Tempo de execução (Jacobi) = %fs\n", timeJacobi);


[xSeidel, erSeidel, iterSeidel, timeSeidel] = fast_seidel(A, b, tol, nMaxIter);
normSeidel = norm(xSeidel, inf);
printf("Norma da solução (Seidel) = %f\n", normSeidel);
printf("Erro (Seidel) = %f\n", erSeidel(iterSeidel));
printf("Número de iterações (Seidel) = %d\n", iterSeidel);
printf("Tempo de execução (Seidel) = %fs\n", timeSeidel);


[xSOR, erSOR, iterSOR, timeSOR] = fast_sor(A, b, tol, nMaxIter, omega);
normSOR = norm(xSOR, inf);
printf("Norma da solução (SOR) = %f\n", normSOR);
printf("Erro (SOR) = %f\n", erSOR(iterSOR));
printf("Número de iterações (SOR) = %d\n", iterSOR);
printf("Tempo de execução (SOR) = %fs\n", timeSOR);

iter = max([iterJacobi, iterSeidel, iterSOR]);

erArrayJacobi = zeros(iter, 1);
for i = 1:iterJacobi
    erArrayJacobi(i) = erJacobi(i);
endfor

erArraySeidel = zeros(iter, 1);
for i = 1:iterSeidel
    erArraySeidel(i) = erSeidel(i);
endfor

erArraySOR = zeros(iter, 1);
for i = 1:iterSOR
    erArraySOR(i) = erSOR(i);
endfor

hf = figure();
plot(1:iter, log(erArrayJacobi), "r", 1:iter, log(erArraySeidel), "g", 1:iter, log(erArraySOR), "b");
xlabel("Número de Iterações");
ylabel("Erro");
legend("Jacobi", "Seidel", "SOR");
title("Número de Iterações X Erro");
print(hf, "rail_1357_iterative.png", "-dpng");
