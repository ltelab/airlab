function val = get_struct_field(S,field)

    if isfield(S,field)
        
        val = getfield(S,field);
        
    else
        
        val = NaN;
        
    end
        
        