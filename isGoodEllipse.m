function isGood = isGoodEllipse(ellipseParam,maxLongAxis,minShortAxis,minArea)
    if (isempty(ellipseParam))
        isGood = false;
    else
        long_axis = ellipseParam.long_axis;
        short_axis = ellipseParam.short_axis;
        isGood = logical(long_axis < maxLongAxis & short_axis > minShortAxis ...
            & pi*long_axis*short_axis > minArea );
    end
end