% Function to sum a series up to n.
% Inputs:
% series - The function handler to the series function.
% x - The x values
% t - The t values
% n - The iteration to sum the series.
% Output:
% S - The total sum of the series.
function S=seriesSum(series,x,t,n)
S=0;
i=1;
    while (i<n+1)
        S=S+series(x,t,i);
        i=i+1;
    end
end