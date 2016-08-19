

x=86;
runs = int64(log(86)/log(3))+1;
coordinates = zeros(runs,1);

for r =1:runs
    disp(strcat('r is ',int2str(r), ' and x is ',int2str(x)));
    ex = (3^(runs-r))
    if(ex>x)
        coordinates(r) = 1
    else
        result = int64(x/ex)
        disp(strcat('adding: ',int2str(result+1)));
        coordinates(r) = result+1;
        x = x - (ex*result)
    end
end
