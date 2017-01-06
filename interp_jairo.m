function dy = interp_jairo(dx,tabx,taby)
% Linear interpolation and extrapolation:
%   dy = interp_jairo(dx,tabx,taby)

dy=NaN(1,length(dx));
    for x = 1:length(dx)
        if isnan(dx(x)), dy(x)=NaN;continue, end
    if dx(x) < tabx(1)
        dy(x) = ((dx(x)-tabx(1))*taby(2)+(tabx(2)-dx(x))*taby(1))/(tabx(2)-tabx(1));
    elseif dx(x) > tabx(end)
        dy(x) = ((dx(x)-tabx(end-1))*taby(end)+(tabx(end)-dx(x))*taby(end-1))/(tabx(end)-tabx(end-1));
    else
        n1 = find(dx(x)>=tabx); n2 = find(dx(x)<=tabx);
        n1 = n1(end); n2 = n2(1);
        if dx(x) == tabx(n1)
            dy(x) = taby(n1);
        else
            dy(x) = ((dx(x)-tabx(n1))*taby(n2)+(tabx(n2)-dx(x))*taby(n1))/(tabx(n2)-tabx(n1));
        end
    end
    end
end