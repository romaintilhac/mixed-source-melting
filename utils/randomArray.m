function [array]=randomArray(iterations,distrib,min_range,max_range,mean,stdev,pd)
%  Generates a random array based on boundaries and distribution parameters

% 0 = uniform (min_range-max_range),
% 1 = normal (based on min_range-max_range),
% 2 = normal (based on the mean and min_range-max_range or mean and stdev)
% 3 = custom (pd must be provided)

if isempty(distrib);clear distrib;end
if isempty(min_range);clear min_range;end
if isempty(max_range);clear max_range;end
if isempty(mean);clear mean;end
if isempty(stdev);clear stdev;end
if isempty(pd);clear pd;end
array=[];

if distrib==0 % uniform (min_range-max_range)
    if exist('min_range','var') && exist('max_range','var')   && ~isnan(min_range)  && ~isnan(max_range)
        array=unifrnd(min_range, max_range,[iterations,1]);
    else
        warning(append('min_range and max_range must be provided')); keyboard;
    end

elseif distrib==1 % normal (based on min_range-max_range)
    if exist('min_range','var') && exist('max_range','var')   && ~isnan(min_range)  && ~isnan(max_range)
        mu = (max_range + min_range) / 2;
        sigma = (max_range - min_range) / 6;
        while numel(array) < iterations
            value = normrnd(mu, sigma);
                if value >= min_range && value <= max_range
                        array = [array; value];
                end
        end
    else
        warning(append('min_range and max_range must be provided for ', input)); keyboard;
    end
    clear mu sigma value

elseif distrib==2 % normal (based on the mean and min_range-max_range or stdev)
    mu = mean;
    if exist('stdev','var') && ~isnan(stdev)
        sigma = stdev;
        if exist('min_range','var') && exist('max_range','var')  && ~isnan(min_range)  && ~isnan(max_range) %to handle both sigma+minmax
            while numel(array) < iterations
                value = normrnd(mu, sigma);
                    if value >= min_range && value <= max_range
                            array = [array; value];
                    end
            end
        end

    elseif exist('min_range','var') && exist('max_range','var')  && ~isnan(min_range)  && ~isnan(max_range)
        sigma = (max_range - min_range) / 6; % stdev
    else
        warning(append('min_range-max_range or stdev must be provided in addition to the mean for ', input)); keyboard;
    end
    while numel(array) < iterations
        value = normrnd(mu, sigma);
        if exist('min_range','var') && exist('max_range','var')  && ~isnan(min_range)  && ~isnan(max_range)
            if value >= min_range && value <= max_range
                    array = [array; value];
            end
        else
            array = [array; value];
        end
    end
    clear mu sigma value

elseif distrib==3 % custom (pd must be provided)
    array=random(pd, iterations, 1);

end

end