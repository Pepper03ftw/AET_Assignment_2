function [LHV_Jkg, LHV_MJkg, stoich] = computeFuelLHVandStoich(fuel, dHf, M)
    % Balance: CxHy + a O2 -> x CO2 + (y/2) H2O
    x = fuel.nC;
    y = fuel.nH;

    a = x + y/4;           % stoich O2 moles per mole fuel (for complete oxidation to CO2 + H2O)
    nCO2 = x;
    nH2O = y/2;

    % Reaction enthalpy (kJ/mol fuel):
    % ΔH = Σ(ν * ΔHf_products) - Σ(ν * ΔHf_reactants)
    dH_rxn_kJmol = (nCO2*dHf.CO2_kJmol + nH2O*dHf.H2Og_kJmol) ...
                   - (fuel.dHf_kJmol + a*dHf.O2_kJmol);

    % LHV is energy released -> positive value
    LHV_kJmol = -dH_rxn_kJmol;

    % Convert to J/kg fuel
    M_fuel = x*M.C + y*M.H;   % kg/mol
    LHV_Jkg = (LHV_kJmol*1000) / M_fuel;
    LHV_MJkg = LHV_Jkg / 1e6;

    % Stoichiometric FAR:
    % For 1 mol fuel: needs a mol O2. Air moles = a / 0.21. Air mass from average molar mass 28.97 g/mol.
    nAir = a / 0.21;
    M_air = 28.97e-3;         % kg/mol
    mAir_per_molFuel = nAir * M_air;
    mFuel_per_molFuel = M_fuel;

    stoich.FAR_stoich = mFuel_per_molFuel / mAir_per_molFuel;
    stoich.O2_mol_per_molFuel = a;
    stoich.air_mol_per_molFuel = nAir;
end