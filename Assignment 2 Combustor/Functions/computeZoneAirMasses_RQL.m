function zone = computeZoneAirMasses_RQL(m3, airSplit)
    % RQL zone air "present in that zone" according to assignment:
    % Rich zone air basis = 0.12*m3 + X1
    % Quench zone air basis = X2
    % Lean/dilution zone air basis = X3 (and/or remaining)
    %
    % If X1/X2/X3 are fractions of m3, then X1 = X1_frac*m3, etc.

    if any(~isfinite([airSplit.rich_frac, airSplit.X1_frac, airSplit.X2_frac, airSplit.X3_frac]))
        warning('X1/X2/X3 fractions not fully specified. Set airSplit.X1_frac, X2_frac, X3_frac.');
    end

    zone.m_air_rich_present   = airSplit.rich_frac*m3 + airSplit.X1_frac*m3;
    zone.m_air_quench_present = airSplit.X2_frac*m3;
    zone.m_air_lean_present   = airSplit.X3_frac*m3;

    zone.m3_total = m3;
end