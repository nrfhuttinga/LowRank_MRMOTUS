function even_integer=make_even(input)
    % Make odd input even by substracting one.
    %
    % Niek Huttinga, UMC Utrecht, 2020.

    input = round(input);

    if mod(input,2)~=0
        even_integer = input - 1 ;
    end
    
end
