function even_integer=make_even(input)
    % Make odd input even by substracting one.
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

    input = round(input);

    if mod(input,2)~=0
        even_integer = input - 1 ;
    else
        even_integer = input;
    end
    
end
