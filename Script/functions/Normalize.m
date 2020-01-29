function new_el = Normalize(whatIsIt, element)

    if(whatIsIt=="matrix")
        new_el=element/element(3,3);
    elseif(whatIsIt=="vector")
        new_el=element/element(3); 
    elseif(whatIsIt=="segment")
        new_el = [Normalize("vector", element(:,1)) ,Normalize("vector", element(:,2))];
    else
        new_el = [];
    end

end

