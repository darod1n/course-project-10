function [P] = progDavlenie(T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    N = length(T);
    P = zeros(1, N);
    for i = 1:N
        if T(1,i) < 1
            P(1,i) = 39 * T(i);
        elseif T(1,i)>= 1 && T(i) < 2
            P(1,i) = 39;
        elseif T(1,i)>= 2 && T(1,i) < 2.5
            y1 = 39; 
            y2 = 29;
            x1 = 2;
            x2 = 2.5;
            k = (y1 - y2) / (x2 - x1);
            b = y1 + k*x1;
            P(1,i) = -k*T(1,i) + b;
        elseif T(1,i)>= 2.5 && T(1,i) < 6
            P(1,i) = 28;
        elseif T(1,i)>= 6 && T(1,i) < 7.5
            y1 = 28; 
            y2 = 100;
            x1 = 6;
            x2 = 7.5;
            k = (y1 - y2) / (x1 - x2);
            b = y1 - k*x1;
            P(1,i) = k*T(1,i) + b;
        elseif T(1,i)>= 7.5 && T(1,i) < 12.5
            P(1,i) = 100;
        elseif T(1,i)>= 12.5 && T(1,i) < 13
            y1 = 100; 
            y2 = 49;
            x1 = 12.5;
            x2 = 13;
            k = (y1 - y2) / (x2 - x1);
            b = y1 + k*x1;
            P(1,i) = -k*T(1,i) + b;
        elseif T(1,i)>= 13 && T(1,i) < 14
            P(1,i) = 49;
        else
            P(i) = 49*14 / T(i);
        end
    end
end

