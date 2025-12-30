function Tad = computeZoneTad_0D(Tin, mf, LHV, cp_gas, zone, m3, airSplit)
    % Compute zone equivalence ratios based on assignment definitions:
    % FAR_rich = mf / (0.12*m3 + X1)
    % FAR_quench = mf / (X2)
    % FAR_lean = mf / (X3)  (if you interpret that way)
    %
    % Then estimate Tad via constant-cp energy balance.
    %
    % NOTE: This is a simplified "available energy / cp" approach and does not explicitly
    % enforce chemical equilibrium or dissociation.

    m_air_rich   = airSplit.rich_frac*m3 + airSplit.X1_frac*m3;
    m_air_quench = airSplit.X2_frac*m3;
    m_air_lean   = airSplit.X3_frac*m3;

    % Total mass flowing through each "zone control volume" (air + fuel),
    % using the mass present there per assignment's zone definition.
    mdot_rich_total   = m_air_rich   + mf;
    mdot_quench_total = m_air_quench + mf;
    mdot_lean_total   = m_air_lean   + mf;

    % Temperature rises (K)
    dT_rich   = (mf * LHV) / (mdot_rich_total   * cp_gas);
    dT_quench = (mf * LHV) / (mdot_quench_total * cp_gas);
    dT_lean   = (mf * LHV) / (mdot_lean_total   * cp_gas);

    Tad.rich_K   = Tin + dT_rich;
    Tad.quench_K = Tin + dT_quench;
    Tad.lean_K   = Tin + dT_lean;
end
