function res = demax(x)
    % Scale input between -1 and 1
    %
    % Tom Bruijnen, UMC Utrecht, 2020.

    res=(x/max(abs(x(:))));

end