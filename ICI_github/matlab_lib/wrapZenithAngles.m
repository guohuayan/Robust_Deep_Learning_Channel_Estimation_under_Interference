function theta = wrapZenithAngles(theta)

    theta = mod(theta,360);
    if (~isempty(coder.target) && isvector(theta))
        % force column indexing as code generation 
        % assumes that vectors are columns
        theta(:) = wrap(theta(:));
    else
        theta = wrap(theta); 
    end
    
    function theta = wrap(theta)
        theta(theta>180) = 360 - theta(theta>180);
    end

end