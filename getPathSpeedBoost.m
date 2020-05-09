function [path, speed, boost] = getPathSpeedBoost(a0, a1, a2, a3, a4, a5, t, N)
%Функция getPathSpeedBoost принимает в качестве аргументов коэффициенты
%a0-a5, t - вектор времени, N - количество точек
    path = zeros(1, N);
    speed = zeros(1, N);
    boost = zeros(1, N);
    for i = 1:N
        path(1, i) = a0 + a1*t(1, i) + a2*t(1, i)^2 + a3*t(1, i)^3 + a4*t(1, i)^4 + a5*t(1, i)^5;
        speed(1, i) = a1 + 2*a2*t(1, i) + 3*a3*t(1, i)^2 + 4*a4*t(1, i)^3 + 5*a5*t(1, i)^4;
        boost(1, i) = 2*a2 + 6*a3*t(1, i) + 12*a4*t(1, i)^2 + 20*a5*t(1, i)^3;
    end
end

