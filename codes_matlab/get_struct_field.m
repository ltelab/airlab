function val = get_struct_field(S,field1,field2,field3)

        if isfield(S,field1)

            if nargin == 2
                val = getfield(S,field1);
                
            elseif nargin == 3 
                val = get_struct_field(getfield(S,field1),field2);
                
            else
                val = get_struct_field(getfield(getfield(S,field1),field2),field3);
                
            end

        else
            val = NaN;

        end   
end
    
    
        
        