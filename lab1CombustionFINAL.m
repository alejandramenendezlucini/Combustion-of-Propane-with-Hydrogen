% Lab 1 Combustion
% Alejandra Menéndez Lucini
% Study on whether adding H2 to a propane based fuel reduces carbon
% emissions

clc
clear all
close all

% Operating conditions
T0 = 298;       % [K]
p0 = 1;         % [atm]
Ru = 8.314;     % [J/mol·K]

% Weights
W_C3H8 = 12*3+8;
W_H2 = 2;
W_O2 = 16*2;
W_N2 = 14*2;
W_CO2 = 12+16*2;
W_CO = 12+16;
W_H2O = 2+16;

% Volumetric hydrogen fraction in fuel mixture
x_H2_pct = [0.0, 0.20, 0.50, 0.80, 1.0];
phi = 0.1:0.1:2.0;

% Carbon (x) and Hydrogen (y) per mole of fuel blend
x = (1 - x_H2_pct) * 3;                     
y = (1 - x_H2_pct) * 8 + x_H2_pct * 2;     
a_st = x + y / 4;             

% NASA-7 POLYNOMIAL COEFFICIENTS
nasa = struct();

% C3H8 (propane)
nasa.C3H8.Tmid = 1000;
nasa.C3H8.high = [ 7.52441520E+00  1.88982820E-02 -6.29210410E-06  9.21614570E-10 -4.86844780E-14 -1.65643940E+04 -1.78383750E+01];
nasa.C3H8.low  = [ 9.28510930E-01  2.64605660E-02  6.03324460E-06 -2.19149530E-08  9.49615440E-12 -1.40579070E+04  1.92255380E+01];

% H2
nasa.H2.Tmid = 1000;
nasa.H2.high = [ 2.99142337E+00  7.00069764E-04 -5.63382869E-08 -9.23157818E-12  1.58275179E-15 -8.35033997E+02 -1.35511017E+00];
nasa.H2.low  = [ 2.34433112E+00  7.98052075E-03 -1.94781510E-05  2.01572094E-08 -7.37611761E-12 -9.17935173E+02  6.83015201E-01];

% O2
nasa.O2.Tmid = 1000;
nasa.O2.high = [ 3.66096065E+00  6.56365811E-04 -1.41149627E-07  2.05797935E-11 -1.29913436E-15 -1.21597718E+03  3.41536279E+00];
nasa.O2.low  = [ 3.78245636E+00 -2.99673416E-03  9.84730201E-06 -9.68129509E-09  3.24372837E-12 -1.06394356E+03  3.65767573E+00];

% N2
nasa.N2.Tmid = 1000;
nasa.N2.high = [ 2.95257637E+00  1.39690040E-03 -4.92631603E-07  7.86010195E-11 -4.60755204E-15 -9.23948688E+02  5.87188762E+00];
nasa.N2.low  = [ 3.53100528E+00 -1.23660988E-04 -5.02999433E-07  2.43530612E-09 -1.40881235E-12 -1.04697628E+03  2.96747038E+00];

% CO2
nasa.CO2.Tmid = 1000;
nasa.CO2.high = [ 4.63651110E+00  2.74145690E-03 -9.95897590E-07  1.60386660E-10 -9.16198570E-15 -4.90249040E+04 -1.93489550E+00];
nasa.CO2.low  = [ 2.35677352E+00  8.98413000E-03 -7.12206320E-06  2.45730080E-09 -1.42885480E-13 -4.83719710E+04  9.90105220E+00];

% CO
nasa.CO.Tmid = 1000;
nasa.CO.high = [ 3.04848590E+00  1.35172810E-03 -4.85794050E-07  7.88536440E-11 -4.69807460E-15 -1.42661170E+04  6.01709770E+00];
nasa.CO.low  = [ 3.57953350E+00 -6.10353690E-04  1.01681430E-06  9.07005860E-10 -9.04424490E-13 -1.43440860E+04  3.50840930E+00];

% H2O
nasa.H2O.Tmid = 1000;
nasa.H2O.high = [ 2.67214600E+00  3.05629300E-03 -8.73026000E-07  1.20099600E-10 -6.39161800E-15 -2.98992100E+04  6.86281700E+00];
nasa.H2O.low  = [ 3.38684200E+00  3.47498200E-03 -6.35469600E-06  6.96858100E-09 -2.50658800E-12 -3.02081100E+04  2.59023300E+00];



% Fuel properties and LHV
H_C3H8 = species_enthalpy('C3H8', T0, nasa, Ru);
H_O2   = species_enthalpy('O2',   T0, nasa, Ru);
H_CO2  = species_enthalpy('CO2',  T0, nasa, Ru);
H_H2O  = species_enthalpy('H2O',  T0, nasa, Ru);

dH_C3H8 = 3*H_CO2 + 4*H_H2O - H_C3H8 - 5*H_O2;   % J/mol_fuel
LHV_C3H8 = -dH_C3H8 / W_C3H8;                     % J/g = kJ/kg
H_H2 = species_enthalpy('H2', T0, nasa, Ru);
dH_H2 = H_H2O - H_H2 - 0.5*H_O2;                  % J/mol_fuel
LHV_H2 = -dH_H2 / W_H2;                           % J/g = kJ/kg

W_fuel   = (1 - x_H2_pct).*W_C3H8 + x_H2_pct.*W_H2;
Y_C3H8   = (1 - x_H2_pct).*W_C3H8 ./ W_fuel;
Y_H2     = x_H2_pct.*W_H2 ./ W_fuel;
LHV_blend = Y_C3H8 * LHV_C3H8 + Y_H2 * LHV_H2;


% Initialize matrices (rows = phi, cols = fuel blend)
T_ad = zeros(length(phi), length(x_H2_pct));
m_CO = zeros(length(phi), length(x_H2_pct));
m_H2_prod = zeros(length(phi), length(x_H2_pct));
m_CO2 = zeros(length(phi), length(x_H2_pct));
m_fuel = zeros(length(phi), length(x_H2_pct));



for i = 1:length(x_H2_pct)

    for j = 1:length(phi)

        H_reac = enthalpy_reactants(T0, phi(j), x_H2_pct(i), x(i), y(i), a_st(i), nasa, Ru);
        res = @(T) enthalpy_products(T, phi(j), x(i), y(i), a_st(i), nasa, Ru, @water_gas_Kp) - H_reac;
        T_ad(j,i) = fzero(res, 2200);
        a = a_st(i) / phi(j);
        Kp = water_gas_Kp(T_ad(j,i), nasa, Ru);
        [b, c, d, e] = rich_composition(a, x(i), y(i), Kp);
        f = a - a_st(i);
        n_total = b + c + d + e + 3.76 * a + f;
        m_prod = b*W_CO2 + c*W_CO + d*W_H2O + e*W_H2 + f*W_O2 + 3.76*a*W_N2;
        m_CO(j,i) = c * W_CO;
        m_H2_prod(j,i) = e * W_H2;
        m_CO2(j,i) = b * W_CO2;
        m_total = m_CO2(j,i) + m_H2_prod(j,i) + m_CO(j,i) + d*W_H2O + f*W_O2 + 3.76*a*W_N2;
        m_CO(j,i) = m_CO(j,i)/m_total;
        m_H2_prod(j,i) = m_H2_prod(j,i)/m_total;
        m_CO2(j,i) = m_CO2(j,i)/m_total;

        m_fuel_1mol = (1-x_H2_pct(i))*W_C3H8 + x_H2_pct(i)*W_H2;
        % mass pollutant per mass fuel
        mCO_per_mfuel(j,i)  = (c * W_CO)  / m_fuel_1mol;
        mCO2_per_mfuel(j,i) = (b * W_CO2) / m_fuel_1mol;

    end

end


spec_CO  = mCO_per_mfuel  ./ (LHV_blend / 1000);   % kg CO  / MJ
spec_CO2 = mCO2_per_mfuel ./ (LHV_blend / 1000);   % kg CO2 / MJ



% ========== PLOTTING ==========
figure('Position', [100, 100, 1050, 420]); % Widescreen aspect ratio for 2 panels
colors = [
    0.00, 0.45, 0.74;  % Pure C3H8 (Blue)
    0.85, 0.33, 0.10;  % 20% H2 (Orange)
    0.93, 0.69, 0.13;  % 50% H2 (Yellow/Gold)
    0.49, 0.18, 0.56;  % 80% H2 (Purple)
    0.47, 0.67, 0.19   % Pure H2 (Green)
];
% Subfigure (a): Adiabatic Flame Temperature
subplot(1, 2, 1)

yyaxis left
colororder(colors);           % Applies custom multi-color cycle to left axis
ax = gca;
ax.YAxis(1).Color = 'k';
p1 = plot(phi, T_ad(:,1), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on
p2 = plot(phi, T_ad(:,2), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
p3 = plot(phi, T_ad(:,3), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
p4 = plot(phi, T_ad(:,4), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
p5 = plot(phi, T_ad(:,5), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
ylabel('T_{ad} [K]')
ylim([500, 2600])


yyaxis right
colororder(colors);          
ax = gca;
ax.YAxis(2).Color = 'k';
p_h2_c3h8 = plot(phi,m_H2_prod(:,1), '-.', 'LineWidth', 1.3);
hold on
p_h2_h2 = plot(phi,m_H2_prod(:,5), '-.', 'LineWidth', 1.3);
p_co_c3h8 = plot(phi, m_CO(:,1), '-.', 'LineWidth',1.3);
p_co_h2 = plot(phi, m_CO(:,5), '-.', 'LineWidth',1.3);
ylabel('Dissociation Mass Fraction')
ylim([0,0.23])


grid on
xlabel('Equivalence ratio \phi')
title('(a) Adiabatic flame temperature')
legend([p1, p2, p3, p4, p5, p_h2_c3h8, p_h2_h2, p_co_c3h8, p_co_h2], ...
       {'Pure C_3H_8 (6-eq)', '0.8C_3H_8 + 0.2H_2', '0.5C_3H_8 + 0.5H_2', ...
        '0.2C_3H_8 + 0.8H_2', 'Pure H_2', 'w_{H_2} (Pure C_3H_8)', 'w_{H_2} (Pure H_2)', 'w_{CO} (Pure C_3H_8)', 'w_{CO} (Pure H_2)'}, ...
       'Location', 'northwest', 'FontSize', 8);
hold off



% Subfigure (b): Specific CO and CO₂ emissions (rich only)
idx_rich = phi > 1;
subplot(1, 2, 2)

% ---- CO (dashed lines) ----
pCO1 = plot(phi(idx_rich), spec_CO(idx_rich,1), '-.', 'Color', colors(1,:), ...
            'LineWidth', 1.5, 'MarkerSize', 5); hold on
pCO2 = plot(phi(idx_rich), spec_CO(idx_rich,2), '-.', 'Color', colors(2,:), ...
            'LineWidth', 1.5, 'MarkerSize', 5);
pCO3 = plot(phi(idx_rich), spec_CO(idx_rich,3), '-.', 'Color', colors(3,:), ...
            'LineWidth', 1.5, 'MarkerSize', 5);
pCO4 = plot(phi(idx_rich), spec_CO(idx_rich,4), '-.', 'Color', colors(4,:), ...
            'LineWidth', 1.5, 'MarkerSize', 5);
pCO5 = plot(phi(idx_rich), spec_CO(idx_rich,5), '-.', 'Color', colors(5,:), ...
            'LineWidth', 1.5, 'MarkerSize', 5);

% ---- CO₂ (solid lines with markers) ----
pCO2_1 = plot(phi(idx_rich), spec_CO2(idx_rich,1), '-o', 'Color', colors(1,:), ...
              'LineWidth', 1.5, 'MarkerSize', 5);
pCO2_2 = plot(phi(idx_rich), spec_CO2(idx_rich,2), '-o', 'Color', colors(2,:), ...
              'LineWidth', 1.5, 'MarkerSize', 5);
pCO2_3 = plot(phi(idx_rich), spec_CO2(idx_rich,3), '-o', 'Color', colors(3,:), ...
              'LineWidth', 1.5, 'MarkerSize', 5);
pCO2_4 = plot(phi(idx_rich), spec_CO2(idx_rich,4), '-o', 'Color', colors(4,:), ...
              'LineWidth', 1.5, 'MarkerSize', 5);
pCO2_5 = plot(phi(idx_rich), spec_CO2(idx_rich,5), '-o', 'Color', colors(5,:), ...
              'LineWidth', 1.5, 'MarkerSize', 5);

grid on
xlabel('Equivalence ratio \phi')
ylabel('Specific emission [kg / MJ]')
title('(b) Specific CO and CO_2 emissions (rich, \phi > 1)')

% Legend – group by species and by blend
legend([pCO1, pCO2, pCO3, pCO4, pCO5, pCO2_1, pCO2_2, pCO2_3, pCO2_4, pCO2_5], ...
       {'CO – Pure C_3H_8', 'CO – 20% H_2', 'CO – 50% H_2', 'CO – 80% H_2', 'CO – Pure H_2', ...
        'CO_2 – Pure C_3H_8', 'CO_2 – 20% H_2', 'CO_2 – 50% H_2', 'CO_2 – 80% H_2', 'CO_2 – Pure H_2'}, ...
       'Location', 'northeast', 'FontSize', 7, 'NumColumns', 2);

hold off



% ========== HELPER FUNCTIONS ==========
function h = species_enthalpy(sp, T, nasa, Ru)
    data = nasa.(sp);
    if T >= data.Tmid
        a = data.high;
    else
        a = data.low;
    end
    hRT = a(1) + a(2)*T/2 + a(3)*T^2/3 + a(4)*T^3/4 + a(5)*T^4/5 + a(6)/T;
    h   = hRT * Ru * T;
end

function H = enthalpy_reactants(T0, phi, x_H2, x, y, a_st, nasa, Ru)
    a = a_st / phi;
    n_C3H8 = 1 - x_H2;
    n_H2   = x_H2;
    n_O2   = a;
    n_N2   = 3.76 * a;
    
    H = n_C3H8 * species_enthalpy('C3H8', T0, nasa, Ru) + ...
        n_H2   * species_enthalpy('H2',   T0, nasa, Ru) + ...
        n_O2   * species_enthalpy('O2',   T0, nasa, Ru) + ...
        n_N2   * species_enthalpy('N2',   T0, nasa, Ru);
end

function H = enthalpy_products(T, phi, x, y, a_st, nasa, Ru, Kp_fun)
    a = a_st / phi;
    if phi <= 1
        b = x;              % CO2
        c = 0;              % CO
        d = y/2;            % H2O
        e = 0;              % H2
        f = a - a_st;       % O2
    else
        Kp = water_gas_Kp(T, nasa, Ru);
        [b, c, d, e] = rich_composition(a, x, y, Kp);
        f = 0;
    end
    n_N2 = 3.76 * a;
    
    H = b * species_enthalpy('CO2', T, nasa, Ru) + ...
        c * species_enthalpy('CO',  T, nasa, Ru) + ...
        d * species_enthalpy('H2O', T, nasa, Ru) + ...
        e * species_enthalpy('H2',  T, nasa, Ru) + ...
        f * species_enthalpy('O2',  T, nasa, Ru) + ...
        n_N2 * species_enthalpy('N2', T, nasa, Ru);
end

function [b, c, d, e] = rich_composition(a, x, y, Kp)
    if x == 0 % Pure H2 case
        b = 0;
        c = 0;
        d = 2 * a;
        e = y/2 - 2*a;
        return;
    end
    fun = @(v) [
        v(1) + v(2) - x;                     % C balance: b + c = x
        v(3) + v(4) - y/2;                   % H balance: d + e = y/2
        2*v(1) + v(2) + v(3) - 2*a;          % O balance: 2b + c + d = 2a
        v(1)*v(4) - Kp*v(2)*v(3)             % Kp = (b*e)/(c*d)
    ];
    v0 = [x*0.6; x*0.4; y/4; y/4];
    options = optimoptions('fsolve', 'Display', 'off');
    sol = fsolve(fun, v0, options);
    b = max(sol(1), 0);
    c = max(sol(2), 0);
    d = max(sol(3), 0);
    e = max(sol(4), 0);
end

function [h, s] = species_thermo(sp, T, nasa, Ru)
    data = nasa.(sp);
    if T >= data.Tmid
        a = data.high;
    else
        a = data.low;
    end
    hRT = a(1) + a(2)*T/2 + a(3)*T^2/3 + a(4)*T^3/4 + a(5)*T^4/5 + a(6)/T;
    h   = hRT * Ru * T;
    sR  = a(1)*log(T) + a(2)*T + a(3)*T^2/2 + a(4)*T^3/3 + a(5)*T^4/4 + a(7);
    s   = sR * Ru;
end

function Kp = water_gas_Kp(T, nasa, Ru)
[H_CO2, S_CO2] = species_thermo('CO2', T, nasa, Ru);
[H_CO,  S_CO]  = species_thermo('CO',  T, nasa, Ru);
[H_H2O, S_H2O] = species_thermo('H2O', T, nasa, Ru);
[H_H2,  S_H2]  = species_thermo('H2',  T, nasa, Ru);  
dG = (H_CO2 - T*S_CO2) + (H_H2 - T*S_H2) - (H_CO - T*S_CO) - (H_H2O - T*S_H2O);
Kp = exp(-dG / (Ru * T));
end