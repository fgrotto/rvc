%% thomas_algorithm.m
%  It calculates a solution for linear systems Ax=d where A is a tridiagonal matrix
function x = thomas_algorithm(A,d)
    [n,m] = size(A);
    a = diag(A,-1);
    b = diag(A);
    c = diag(A,1);
    
    assert(n==m, "The A matrix must be nxn")
    assert(n==length(d), "The dimention of d must be n with A nxn")
  
    % forward elimination
    for k=2:1:n 
      m = a(k-1)/b(k-1);
      b(k) = b(k)-m*c(k-1);
      d(k) = d(k)-m*d(k-1);
    end

    % backward substitution
    x(n) = d(n)/b(n);
    for k = n-1:-1:1
      x(k) = (d(k)-c(k)*x(k+1))/b(k);
    end
end
