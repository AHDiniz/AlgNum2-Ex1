load "685_bus.mat"

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

[xJacobi, erJacobi, iterJacobi] = fast_jacobi(A, b, tol, nMaxIter);
normJacobi = norm(xJacobi, inf);
printf("Norma da solução (Jacobi) = %f\n", normJacobi);
printf("Erro (Jacobi) = %f\n", erJacobi(iterJacobi));
printf("Número de iterações (Jacobi) = %d\n", iterJacobi);

[xSeidel, erSeidel, iterSeidel] = fast_seidel(A, b, tol, nMaxIter);
normSeidel = norm(xSeidel, inf);
printf("Norma da solução (Seidel) = %f\n", normSeidel);
printf("Erro (Seidel) = %f\n", erSeidel(iterSeidel));
printf("Número de iterações (Seidel) = %d\n", iterSeidel);

[xSOR, erSOR, iterSOR] = fast_sor(A, b, tol, nMaxIter, omega);
normSOR = norm(xSOR, inf);
printf("Norma da solução (SOR) = %f\n", normSOR);
printf("Erro (SOR) = %f\n", erSOR(iterSOR));
printf("Número de iterações (SOR) = %d\n", iterSOR);
