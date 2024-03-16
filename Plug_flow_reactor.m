function PFR_profiles
    % Given data
    F_A0 = 10; % mol/s
    V0 = 10; % dm^3/s
    V = 1000; % dm^3
    T0 = 350; % K
    Ua = 10; % cal/s.dm^3.K
    Ta = 373.15; % K (example ambient temperature in K, assuming 100C)
    % Heat capacities (cal/mol.K)
    CpA = 30; CpB = 30; CpC = 30; CpD = 90;
    % Heat of reactions (cal/mol)
    dHr1 = -50000; dHr2 = 5000;
    % Rate constants
    k1 = 0.043; % dm^3/mol.s
    k2 = 0.4; % dm^3/mol.s

    % Initial conditions
    F_B0 = F_A0; % Equimolar in A and B
    F_C0 = 0;
    F_D0 = 0;
    y0 = [F_A0, F_B0, F_C0, F_D0, T0];

    % Volume span
    Vspan = [0 V];

    % Solve ODEs
    [V, Y] = ode45(@(V, y) reactorODE(V, y, k1, k2, Ua, Ta, CpA, CpB, CpC, CpD, dHr1, dHr2), Vspan, y0);

    % Extract solutions
    F_A = Y(:, 1);
    F_B = Y(:, 2);
    F_C = Y(:, 3);
    F_D = Y(:, 4);
    T = Y(:, 5);

    % Plot results
    figure;
    subplot(2,2,1);
    plot(V, T);
    title('Temperature Profile');
    xlabel('Volume (dm^3)');
    ylabel('Temperature (K)');

    subplot(2,2,2);
    plot(V, [F_A, F_B, F_C, F_D]);
    title('Flow Rate Profiles');
    xlabel('Volume (dm^3)');
    ylabel('Flow Rate (mol/s)');
    legend('F_A', 'F_B', 'F_C', 'F_D');

    % Assuming constant volumetric flow rate for concentration profiles
    C_A = F_A / V0;
    C_B = F_B / V0;
    C_C = F_C / V0;
    C_D = F_D / V0;

    subplot(2,2,3);
    plot(V, [C_A, C_B, C_C, C_D]);
    title('Concentration Profiles');
    xlabel('Volume (dm^3)');
    ylabel('Concentration (mol/dm^3)');
    legend('C_A', 'C_B', 'C_C', 'C_D');
end

function dydV = reactorODE(V, y, k1, k2, Ua, Ta, CpA, CpB, CpC, CpD, dHr1, dHr2)
    F_A = y(1);
    F_B = y(2);
    F_C = y(3);
    F_D = y(4);
    T = y(5);

    % Reaction rates
    r1 = k1 * (F_A/V0) * (F_B/V0);
    r2 = k2 * (F_B/V0) * (F_C/V0);

    % Differential equations
    dFAdV = -r1;
    dFBdV = -r1 - 2*r2;
    dFCdV = 2*r1 - r2;
    dFDdV = r2;
    dTdV = (-(dHr1 * r1 + dHr2 * r2) + Ua * (Ta - T)) / (F_A*CpA + F_B*CpB + F_C*CpC + F_D*CpD);

    dydV = [dFAdV; dFBdV; dFCdV; dFDdV; dTdV];
end
