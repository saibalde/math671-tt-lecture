function [ttf, XEVALS, FEVALS] = tt_fun(fun, relAcc, d, type)
%TT_FUN Construct TT approximation to function using TT
%
%   [ttf, xevals, fevals] = tt_fun(fun, relAcc, d, type) computes the TT
%   approximation of the vector generated by evauating function handle fun
%   at centers of the cells generated by dividing interval [0, 1] into
%   2^d equal segments. The TT decomposition ttf is then obtained by
%   running
%
%       *   TT-SVD if type is 'svd'
%       *   TT-AMEn-cross if type is 'amen-cross'
%
%   with relative accuracy relAcc on the tensorized function value vector.
%   xevals and fevals record where the function was evaluated.

    h = 1 / 2^d;
    x = (0.5 * h : h : 1.0)';

    if strcmp(type, 'svd')
        funvals = fun(x);
        funvals = reshape(funvals, 2 * ones(1, d));
        ttf = tt_tensor(funvals, relAcc);

        XEVALS = x;
        FEVALS = funvals;
    elseif strcmp(type, 'amen-cross')
        kickRank = 4;
        maxRank  = 400;
        numSweep = ceil(maxRank / kickRank);

        XEVALS = [];
        FEVALS = [];

        ttf = amen_cross(2 * ones(1, d), @fun_eval, relAcc, ...
                         'nswp', numSweep, 'kickrank', kickRank, ...
                         'tol_exit',relAcc);
        ttf = round(ttf, 1.0e-15);
    else
        error('type must be either svd or amen-cross');
    end

    function vals = fun_eval(inds)
        n = size(inds, 1);
        vals = zeros(n);

        for k = 1 : n
            i = 1;
            for l = 1 : d
                i = i + 2^(l - 1) * (inds(k, l) - 1);
            end

            vals(k) = fun(x(i));

            XEVALS = [XEVALS; x(i)];
            FEVALS = [FEVALS; vals(k)];
        end
    end
end