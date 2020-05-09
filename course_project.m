clc;
clear variables;
close all;

%% Начальные условия
T = 26;
g = -9.8;
dt = 0.01;
t = 0:dt:T;
N = numel(t);
C = [0 1];
D = 0;
B = [0; 1 ];
l1=-8;
l2=-15;

A=[
    0 1;
    0 0
];

% sys = ss(A, B, C, D)
% sys2 = mksys(A, B, C, D);

%% Задание полюсов 
% p = [-8; -15];
% pe = 2*[-8; -7];
% K = place(A, B, p);
% L = acker(A, B, pe)';

%% Формирование наблюдателя с выбранными полюсами
% newsys = append(sys, K)
% rsys = reg(sys, K, L);
% sys1 = parallel(sys, rsys, [1], [1], [2]);
% sys2 = connect(sys1, );

%% Проверка управляемости и наблюдаемости
% Co = ctrb(A, B);
% Unco = length(A) - rank(Co);
% Do = obsv(sys);
% Undo = length(A) - rank(Do);

%%  Расчёт траектории по эта, скорости и ускорения с учётом алгоритмов описывающих функций

%   Начальные условия по координате Кси
xiPathInit = 100;
xiSpeedInit = 0;
xiBoostInit = 0;

%   Конечные условия по координате Кси
xiPathFin = 0;
xiSpeedFin = 0;
xiBoostFin = 0;

%   Получение коэффициентов по координате Кси
[coefXi0, coefXi1, coefXi2, coefXi3, coefXi4, coefXi5] = getCoefficient(xiPathInit, xiSpeedInit, xiBoostInit, xiPathFin, xiSpeedFin, xiBoostFin, T);

%   pathXi - траектория пути по Кси
%   speedXi - скорость по Кси
%   boostXi - ускорения по Кси
[pathXi, speedXi, boostXi] = getPathSpeedBoost(coefXi0, coefXi1, coefXi2, coefXi3, coefXi4, coefXi5, t, N);

%   Начальные условия по координате Эта
etaPathInit = 800;
etaSpeedInit = -70;
etaBoostInit = etaSpeedInit ^ 2 / (2 * etaPathInit);

%   Конечные условия по координате Эта
etaPathFin = 20;
etaSpeedFin = -2;
etaBoostFin = 0;

%   Получение коэффициентов по координате Эта
[coefEta0, coefEta1, coefEta2, coefEta3, coefEta4, coefEta5] = getCoefficient(etaPathInit, etaSpeedInit, etaBoostInit, etaPathFin, etaSpeedFin, etaBoostFin, T);

%   pathEta - траектория пути по Эта
%   speedEta - скорость по Эта
%   boostEta - ускорения по эта
[pathEta, speedEta, boostEta] = getPathSpeedBoost(coefEta0, coefEta1, coefEta2, coefEta3, coefEta4, coefEta5, t, N);

%   Начальные условия по координат Дзета
zetaPathInit = -100;
zetaSpeedInit = 0;
zetaBoostInit = 0;

%   Конечные условия по координате Дзета
zetaPathFin = 0;
zetaSpeedFin = 0;
zetaBoostFin = 0;

%   Получение коэффициентов по координате Дзета
[coefZeta0, coefZeta1, coefZeta2, coefZeta3, coefZeta4, coefZeta5] = getCoefficient(zetaPathInit, zetaSpeedInit, zetaBoostInit, zetaPathFin, zetaSpeedFin, zetaBoostFin, T);

%   pathZeta - траектория пути по Кси
%   speedZeta - скорость по Кси
%   boostZeta - ускорения по Кси
[pathZeta, speedZeta, boostZeta] = getPathSpeedBoost(coefZeta0, coefZeta1, coefZeta2, coefZeta3, coefZeta4, coefZeta5, t, N);

%%  Моделирование системы по желаемой траектории
% f0 = [y01 y0]
% [f, t, h] = lsim(sys, y, t, f0);
% [freg, t, hreg ] = lsim(rsys, y, t, f0);

%%  Графики траекторий
%   Траектория по Кси
figure('Name', 'Траектория движения по координате Кси');
plot(t, pathXi);
title('Траектория движения по координате Кси');

%   Траектория по Эта
figure('Name', 'Траектория движения по координате Эта');
plot(t, pathEta);
title('Траектория движения по координате Эта')

%   Траектория по Дзета
figure('Name', 'Траектория движения по координате Дзета');
plot(t, pathZeta);
title('Траектория движения по координате Дзета')

%   Трехмерный график
figure('Name', 'Траектория движения в 3D');
plot3(pathZeta, pathXi, pathEta);

%%  Графики изменения скорости
%   Скорость по Кси
figure('Name', 'График изменения скорости')
plot(t, speedXi);
title('График изменения скорости по Кси');

%   Скорость по Эта
figure('Name', 'График изменения скорости')
plot(t, speedEta);
title('График изменения скорости по Эта');

%   Скорость по Дзета
figure('Name', 'График изменения скорости')
plot(t, speedZeta);
title('График изменения скорости по Дзета');

%%  Графики изменения ускорения
%   Ускорение по Кси
figure('Name', 'График изменения ускорения');
plot(t, boostXi);
title('График изменения ускорения по Кси');

%   Ускорение по Эта
figure('Name', 'График изменения ускорения');
plot(t, boostEta);
title('График изменения ускорения по Эта');

%   Ускорение по Дзета
figure('Name', 'График изменения ускорения');
plot(t, boostZeta);
title('График изменения ускорения по Дзета');

%%
%{
%   Реакция системы на импульс
figure('Name', 'Импульс');
grid on;
impulse(sys);

%   Диаграмма Боде
figure('Name', 'Диаграмма Боде');
bode(sys)

%   Реакция системы на единицу
figure('Name', 'Реакция системы на единицу')
step(sys)

%   График траектор
figure('Name', 'Траектория')
plot(t, y, t, f)
%}