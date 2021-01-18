function SkipLegendEntry(objects)
    
    for i=1:numel(objects)
        objects{i}.HandleVisibility = 'off';
    end
    
end

