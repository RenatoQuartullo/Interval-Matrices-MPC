function P = build_permutation_matrix(N, n, m)
    % N: number of x(i) vectors is N+1
    %    number of u(i) vectors is N
    % n: dimension of x(i)
    % m: dimension of u(i)
    %
    % Y = [x(0); ...; x(N); u(0); ...; u(N-1)]
    % PY = [x(0); u(0); x(1); u(1); ...; x(N-1); u(N-1); x(N)]

    total_input_dim = (N+1)*n + N*m;
    total_output_dim = total_input_dim;

    P = zeros(total_output_dim);

    for k = 0:N-1
        x_out_start = k*(n + m) + 1;
        u_out_start = x_out_start + n;

        x_in_start = k*n + 1;
        u_in_start = (N+1)*n + k*m + 1;

        P(x_out_start:x_out_start+n-1, x_in_start:x_in_start+n-1) = eye(n);
        P(u_out_start:u_out_start+m-1, u_in_start:u_in_start+m-1) = eye(m);
    end

    % Last x(N) goes at the end
    xN_out_start = N*(n + m) + 1;
    xN_in_start = N*n + 1;
    P(xN_out_start:xN_out_start+n-1, xN_in_start:xN_in_start+n-1) = eye(n);
end
