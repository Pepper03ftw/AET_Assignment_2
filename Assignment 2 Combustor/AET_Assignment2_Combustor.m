%% AE4238 Assignment 2 - RQL Combustor (0-D framework)
% You must paste in P3, T3, T4, mdot3 values for the thrust points from the provided plots.
% All calculations follow the assignment statement.

clear; clc; close;

%% ---------------------------
% 0) USER INPUTS (FROM PLOTS)
% ---------------------------
% Thrust settings [kN]
thrust_kN = [0, 50, 100, 200, 300];  % "0" used for idle; replace if you use a different idle thrust reference

% Paste values from plots (same length as thrust_kN)
% P3: HPC exit pressure [kPa] (from "Variation of P3" plot)
P3_kPa  = [NaN, NaN, NaN, NaN, NaN];

% T3: HPC exit temperature [K] (from "Variation of T3" plot)
T3_K    = [NaN, NaN, NaN, NaN, NaN];

% T4: Burner exit temperature [K] (from "Variation of T4" plot)
T4_K    = [NaN, NaN, NaN, NaN, NaN];

% mdot3: HPC exit flow W3 [kg/s] (from "Mass Flow Variation" plot)
mdot3_kgps = [NaN, NaN, NaN, NaN, NaN];

% Quick sanity check:
assert(numel(thrust_kN)==numel(P3_kPa) && numel(P3_kPa)==numel(T3_K) && ...
       numel(T3_K)==numel(T4_K) && numel(T4_K)==numel(mdot3_kgps), ...
       'All input arrays must have the same length.');

%% ---------------------------
% 1) GIVEN CONSTANTS / DATA
% ---------------------------
% Specific heats (given)
cp_air = 1000;    % J/kgK
cp_gas = 1150;    % J/kgK

% Gamma values (not used directly in this 0-D energy balance, kept for completeness)
gamma_air = 1.4;
gamma_gas = 1.33;

% Universal / gas constants
R_u = 8.314462618;             % J/mol/K
M_air = 28.97e-3;              % kg/mol (approx.)
R_air = R_u / M_air;           % J/kg/K

% Fuel assumption: SAF as C11H23 (assignment)
fuel.nC = 11;
fuel.nH = 23;
fuel.dHf_kJmol = -340;          % kJ/mol (given)

% Enthalpies of formation (kJ/mol), given
% Use gas-phase H2O as requested/typical for LHV
dHf.CO2_kJmol = -393.0;
dHf.H2Og_kJmol = -242.0;
dHf.CO_kJmol  = -110.5;
dHf.O2_kJmol  = 0.0;
dHf.N2_kJmol  = 0.0;

% Molar masses (kg/mol)
M.C = 12.00e-3;
M.H = 1.008e-3;
M.O = 16.00e-3;

% Residence time (given)
tau = 6e-3; % s

% RQL staging assumptions from assignment slides:
airSplit.rich_frac = 0.12;   % 12% of m3 goes to rich zone (given on schematic)
% Additional splits X1, X2, X3 depend on your combustor schematic interpretation.
% If your brief expects specific numbers, define them here.
airSplit.X1_frac = NaN;  % fraction of m3 that enters in addition to 0.12*m3 in rich zone
airSplit.X2_frac = NaN;  % fraction of m3 that enters at quench
airSplit.X3_frac = NaN;  % fraction of m3 that enters at dilution/lean

% NOTE: If X1/X2/X3 are not provided numerically, you must make assumptions and justify them.

%% ---------------------------
% 2) FUEL LHV (CALORIFIC VALUE)
% ---------------------------
[LHV_Jkg, LHV_MJkg, stoich] = computeFuelLHVandStoich(fuel, dHf, M);

fprintf('Fuel LHV (based on H2O(g)) = %.3f MJ/kg\n', LHV_MJkg);
fprintf('Stoichiometric FAR = %.6f (kg fuel / kg air)\n', stoich.FAR_stoich);

%% ---------------------------
% 3) FUEL FLOW AT EACH THRUST
% ---------------------------
% Simple combustor energy balance:
% mdot_f * LHV = mdot_air * cp_gas * (T4 - T3)
% (neglecting fuel sensible enthalpy, pressure losses, dissociation, etc.)
dT = T4_K - T3_K;
mdot_f_kgps = (mdot3_kgps .* cp_gas .* dT) ./ LHV_Jkg;

% Optional: guard against negative/NaN due to missing inputs
mdot_f_kgps(~isfinite(mdot_f_kgps)) = NaN;

%% ---------------------------
% 4) OVERALL EQUIVALENCE RATIO
% ---------------------------
FAR_actual = mdot_f_kgps ./ mdot3_kgps;
phi_overall = FAR_actual ./ stoich.FAR_stoich;

%% ---------------------------
% 5) ZONE AIR AMOUNTS @ 300 kN
% ---------------------------
idx300 = find(thrust_kN == 300, 1);
if isempty(idx300)
    warning('300 kN not found in thrust_kN array; skipping zone split calculations.');
    zone = [];
else
    m3_300 = mdot3_kgps(idx300);
    zone = computeZoneAirMasses_RQL(m3_300, airSplit);
end

%% ---------------------------
% 6) ADIABATIC FLAME TEMPERATURES (RICH / QUENCH / LEAN)
% ---------------------------
% This is a simplified 0-D constant-cp estimate:
% Tad = Tin + (mdot_f * LHV) / (mdot_zone_total * cp_gas)
% where "mdot_zone_total" includes the mass flow present in that zone (air + upstream products).
% To respect the assignment definition of zone FAR:
% FAR_rich   = mdot_f / (0.12*m3 + X1)
% FAR_quench = mdot_f / (X2)
% Lean zone similarly uses remaining air added in dilution.

if ~isempty(zone) && isfinite(mdot_f_kgps(idx300))
    Tin = T3_K(idx300);        % assume T3 is inlet to combustor/rich zone
    mf  = mdot_f_kgps(idx300);

    Tad = computeZoneTad_0D(Tin, mf, LHV_Jkg, cp_gas, zone, m3_300, airSplit);
else
    Tad = [];
end

%% ---------------------------
% 7) COMBUSTOR VOLUME FROM RESIDENCE TIME
% ---------------------------
% V = (mdot / rho) * tau
% rho at combustor inlet approximated using P3, T3 (ideal gas)
rho3 = (P3_kPa*1000) ./ (R_air .* T3_K);     % kg/m^3
V_m3 = (mdot3_kgps ./ rho3) * tau;
V_m3(~isfinite(V_m3)) = NaN;

%% ---------------------------
% 8) HEAT DENSITY AT EACH THRUST
% ---------------------------
% Heat release rate: Qdot = mdot_f * LHV [W]
Qdot_W = mdot_f_kgps .* LHV_Jkg;

% Heat density (MW/m^3)
HD_MW_m3 = (Qdot_W ./ V_m3) / 1e6;

% Heat density normalized by pressure (MW/m^3/bar)
P3_bar = P3_kPa / 100; % since 1 bar = 100 kPa
HD_MW_m3_bar = HD_MW_m3 ./ P3_bar;

%% ---------------------------
% 9) OUTPUT TABLE
% ---------------------------
Results = table(thrust_kN(:), P3_kPa(:), T3_K(:), T4_K(:), mdot3_kgps(:), ...
                mdot_f_kgps(:), FAR_actual(:), phi_overall(:), V_m3(:), ...
                HD_MW_m3(:), HD_MW_m3_bar(:), ...
    'VariableNames', {'Thrust_kN','P3_kPa','T3_K','T4_K','mdot3_kgps', ...
                      'mdot_f_kgps','FAR_actual','phi_overall','V_m3', ...
                      'HeatDensity_MW_m3','HeatDensity_MW_m3_bar'});

disp(Results);

% Optional plots (uncomment if useful)
% figure; plot(thrust_kN, mdot_f_kgps, 'o-'); xlabel('Thrust [kN]'); ylabel('Fuel flow [kg/s]'); grid on;
% figure; plot(thrust_kN, phi_overall, 'o-'); xlabel('Thrust [kN]'); ylabel('\phi overall [-]'); grid on;
% figure; plot(thrust_kN, HD_MW_m3, 'o-'); xlabel('Thrust [kN]'); ylabel('Heat Density [MW/m^3]'); grid on;
