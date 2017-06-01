function [] = show_tomo(A,idx,timedelay)

is_matrix = isnumeric(A);

if nargin < 2 || isempty(idx)
    if is_matrix
        idx = 1:size(A,1);
    else
        idx = A([],'size');
        idx = 1:idx(1);
    end
end

if nargin < 3 || isempty(timedelay)
    timedelay = 0.1;
end

if is_matrix
    ca = full([min(A(:)),max(A(:))]);
else
    ca = [0,eps];
end

if is_matrix
    M = size(A,1);
    N = sqrt(size(A,2));
else
    MN = A([],'size');
    M = MN(1);
    N = sqrt(MN(2));
end

h = imagesc(zeros(N,N));
axis image
colorbar
caxis(ca)

if is_matrix
    for k = idx
        set(h,'CData',reshape(A(k,:),N,N))
        title(sprintf('Row %d of %d',k,M));
        pause(timedelay)
    end
else
    for k = idx
        e = zeros(M,1);
        e(k) = 1;
        ray = A(e,'transp');
        ca = [min([ca(1);ray]),max([ca(2);ray])];
        set(h,'CData',reshape(ray,N,N));
        title(sprintf('Row %d of %d',k,M));
        caxis(ca);
        pause(timedelay)
    end
end