function x = back_sub(A, b)
    % Solves the system Ax = b if A is upper triangular
    if size(A,1) ~= size(A,2)
        error("Input matrix A must be square.")
    end
    if any(any(triu(A) ~= A))
        error("Input matrix A must be upper triangular");
    end
    x = zeros(size(b));

    for i = length(b):-1:1
        rhs = b(i);
        for j = i+1:length(b)
            rhs = rhs - A(i,j)*x(j);
        end
        x(i) = rhs/A(i,i);
    end
end