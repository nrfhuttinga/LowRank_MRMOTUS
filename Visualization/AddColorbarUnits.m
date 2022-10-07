function AddColorbarUnits( colobar_object , units)
% set(colobar_object,'FontSize',colobar_object.FontSize-2);
    for i=1:length(colobar_object.TickLabels)
        colobar_object.TickLabels{i}=[colobar_object.TickLabels{i},' ',units];
    end

end