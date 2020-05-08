function [coef0, coef1, coef2, coef3, coef4, coef5] = getCoefficient(pathInit, speedInit, boostInit, pathFin, speedFin, boostFin, T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    coef0 = pathInit;
    coef1 = speedInit;
    coef2 = boostInit/2;
    coef3 = boostFin /(2 * T) - (6 * speedInit) / T ^ 2 - (3 * boostInit)/(2*T) - (4*speedFin)/T^2 - (10*pathInit)/T^3 + (10*pathFin)/T^3;
    coef4 = (15 * pathInit)/T^4 + (8 * speedInit) / T ^ 3 + (3 * boostInit)/(2*T^2) + (7*speedFin)/T^3 - boostFin/T^2 - (15*pathFin)/T ^ 4;
    coef5 = boostFin / (2 * T ^ 3) - (3 * speedInit) / T ^ 4 - boostInit / (2 * T ^ 3) - (3 * speedFin) / T ^ 4 - (6 * pathInit) / T ^ 5 + (6 * pathFin) / T ^ 5;
end

