function result = tensorFroNorm(ten)
    result = 0;
    m = size(ten,1);
    n = size(ten,2);
    z = size(ten,3);
    for i=1:m
        for j=1:n
            for k=1:z
                result = result + abs(ten(i,j,k)^2);
            end
        end
    end
    result = sqrt(result);
end