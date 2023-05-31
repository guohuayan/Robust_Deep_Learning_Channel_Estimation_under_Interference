function [idx] = pilot_idx(num_pilot)
    if num_pilot==1
        idx=3;
    elseif num_pilot==2
        idx=[3,12];
    elseif num_pilot==3
        idx=[3,8,12];
    elseif num_pilot==4
        idx=[3,6,9,12];
    end
end

