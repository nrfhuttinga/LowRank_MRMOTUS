function varargout = mydeal(x)
for k = 1:length(x)
    varargout{k} = x{k};
end
end