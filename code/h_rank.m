function R = h_rank(K, L)
%H_RANK Compute the hierarchical structure of ranks
%
%   h_rank(K, L) computes and returns the hierarchical rank structure of
%   input matrix K at level L depth of the hierarchy. K must be a 2^d x 2^d
%   matrix with integer d, and d >= L >= 2.

N = size(K, 1);

if mod(N, 2^L) ~= 0
    error('N and L values lead to non-uniform hierarchy');
end

R = diagonalBlock(K, L);

figure(1);
clf;
image(log10(K), 'CDataMapping', 'scaled');
title('log_{10}(K_{ij})', 'fontsize', 18);
ax = gca;
ax.FontSize = 16;
axis square;
colorbar;

figure(2);
clf;
image(R, 'CDataMapping', 'scaled');
title('h-rank(K_{ij}), level = ' + string(L), 'fontsize', 18);
ax = gca;
ax.FontSize = 16;
axis square;
colorbar;

end

function R = diagonalBlock(K, L)

n = size(K, 1);
R = zeros(n, n);

if L < 2
    error('Level should not be at least 2');
elseif L == 2
    start = 1 + (0 : 3) * n / 4;
    stop = (1 : 4) * n / 4;
    
    for i = 1 : 4
        rowInds = start(i) : stop(i);
        for j = 1 : 4
            colInds = start(j) : stop(j);
            R(rowInds, colInds) = rank(K(rowInds, colInds));
        end
    end
else
    start = 1 + (0 : 1) * n / 2;
    stop = (1 : 2) * n / 2;
    
    for i = 1 : 2
        rowInds = start(i) : stop(i);
        for j = 1 : 2
            colInds = start(j) : stop(j);
            if i == j
                R(rowInds, colInds) = ...
                    diagonalBlock(K(rowInds, colInds), L - 1);
            elseif i > j
                R(rowInds, colInds) = ...
                    lowerDiagonalBlock(K(rowInds, colInds), L - 1);
            elseif i < j
                R(rowInds, colInds) = ...
                    upperDiagonalBlock(K(rowInds, colInds), L - 1);
            end
        end
    end
end

end

function R = upperDiagonalBlock(K, L)

n = size(K, 1);
R = zeros(n, n);

if L < 2
    error('Level should not be at least 2');
elseif L == 2
    start = 1 + (0 : 3) * n / 4;
    stop = (1 : 4) * n / 4;
    
    for i = 3 : 4
        rowInds = start(i) : stop(i);
        for j = 1 : 2
            colInds = start(j) : stop(j);
            R(rowInds, colInds) = rank(K(rowInds, colInds));
        end
    end
    
    rowInds = start(1) : stop(2);
    colInds = start(1) : stop(2);
    R(rowInds, colInds) = rank(K(rowInds, colInds));
    
    rowInds = start(1) : stop(2);
    colInds = start(3) : stop(4);
    R(rowInds, colInds) = rank(K(rowInds, colInds));
    
    rowInds = start(3) : stop(4);
    colInds = start(3) : stop(4);
    R(rowInds, colInds) = rank(K(rowInds, colInds));
else
    start = 1 + (0 : 1) * n / 2;
    stop = (1 : 2) * n / 2;
    
    for i = 1 : 2
        rowInds = start(i) : stop(i);
        for j = 1 : 2
            colInds = start(j) : stop(j);
            if i == 2 && j == 1
                R(rowInds, colInds) = ...
                    upperDiagonalBlock(K(rowInds, colInds), L - 1);
            else
                R(rowInds, colInds) = rank(K(rowInds, colInds));
            end
        end
    end
end

end

function R = lowerDiagonalBlock(K, L)

n = size(K, 1);
R = zeros(n, n);

if L < 2
    error('Level should not be at least 2');
elseif L == 2
    start = 1 + (0 : 3) * n / 4;
    stop = (1 : 4) * n / 4;
    
    for i = 1 : 2
        rowInds = start(i) : stop(i);
        for j = 3 : 4
            colInds = start(j) : stop(j);
            R(rowInds, colInds) = rank(K(rowInds, colInds));
        end
    end
    
    rowInds = start(1) : stop(2);
    colInds = start(1) : stop(2);
    R(rowInds, colInds) = rank(K(rowInds, colInds));
    
    rowInds = start(3) : stop(4);
    colInds = start(1) : stop(2);
    R(rowInds, colInds) = rank(K(rowInds, colInds));
    
    rowInds = start(3) : stop(4);
    colInds = start(3) : stop(4);
    R(rowInds, colInds) = rank(K(rowInds, colInds));
else
    start = 1 + (0 : 1) * n / 2;
    stop = (1 : 2) * n / 2;
    
    for i = 1 : 2
        rowInds = start(i) : stop(i);
        for j = 1 : 2
            colInds = start(j) : stop(j);
            if i == 1 && j == 2
                R(rowInds, colInds) = ...
                    lowerDiagonalBlock(K(rowInds, colInds), L - 1);
            else
                R(rowInds, colInds) = rank(K(rowInds, colInds));
            end
        end 
    end
end

end