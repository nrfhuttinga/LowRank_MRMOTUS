function res = c12d(vec)
    % Fill dimensions vectors to 12D
    % Tom Bruijnen, UMC Utrecht, 2020
    
    res=vec;
    res(end+1:12)=1;

    % END
end