function y = Afun_matrix(x,transp_flag,A)
switch transp_flag
    case 'size'
        y = size(A);
    case 'notransp'
        y = A*x;
    case 'transp'
        y = A'*x;
    otherwise
        error('transp_flag must be ''size'', ''notransp'' or ''transp''')
end
