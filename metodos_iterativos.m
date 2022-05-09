load "685_bus.mat"

A = Problem.A;
n = rows(A);

solution = ones(n, 1);
b = A * solution;

isDiagDom = diagonal_dominante(A);
printf("É diagonal dominante? %s\n", isDiagDom);

tol = 0.000001;
printf("Tolerância = %f\n", tol);
nMaxIter = 2000;
printf("Número máximo de iterações = %f\n", nMaxIter);
omega = 1.5;
printf("Omega = %f\n", omega);

[xJacobi, erJacobi, iterJacobi] = jacobi(A, b, tol, nMaxIter);
normJacobi = norm(xJacobi, 2);
printf("Norma da solução (Jacobi) = %f\n", normJacobi);
printf("Erro (Jacobi) = %f\n", erJacobi);
printf("Número de iterações (Jacobi) = %f\n", iterJacobi);

[xSeidel, erSeidel, iterSeidel] = seidel(A, b, tol, nMaxIter);
normJacobi = norm(xSeidel, 2);
printf("Norma da solução (Seidel) = %f\n", normSeidel);
printf("Erro (Seidel) = %f\n", erSeidel);
printf("Número de iterações (Seidel) = %f\n", iterSeidel);

[xSOR, erSOR, iterSOR] = sor(A, b, tol, nMaxIter, omega);
normJacobi = norm(xSOR, 2);
printf("Norma da solução (SOR) = %f\n", normSOR);
printf("Erro (SOR) = %f\n", erSOR);
printf("Número de iterações (SOR) = %f\n", iterSOR);
