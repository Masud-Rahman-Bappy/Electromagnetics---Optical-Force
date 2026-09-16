function field = BesselBeamGasLiquidForceFunciton(beam, grid, zM, geometryMode)
arguments
    beam (1,1) struct
    grid (1,1) struct
    zM (1,1) double {mustBeFinite, mustBeReal}
    geometryMode (1,1) string {mustBeMember(geometryMode,["legacy","consistent"])}
end

%% 1. Coordinates -- replace the original nested x/y loops with array algebra.
X = grid.X;
Y = grid.Y;
transverseX = zM*sin(beam.tiltRad) + X*cos(beam.tiltRad);
rho = hypot(transverseX, Y);
axialM = zM*cos(beam.tiltRad) - X*sin(beam.tiltRad);
if geometryMode == "legacy"
    % Preserve the source azimuth convention and transpose on square grids.
    phi = atan(Y ./ transverseX) + pi*(X < 0) + 2*pi*(X > 0 & Y < 0);
    center = X == 0 & Y == 0;
    phi(center) = 0;
    rho(center) = 0;
    axialM(center) = 0;
    phi = phi.';
else
    % atan2 treats all quadrants and the axis without manual branch logic.
    phi = atan2(Y, transverseX);
end

%% 2. Bessel terms -- use analytic limits at rho=0 instead of a tiny divisor.
m = beam.order;
u = beam.q*rho;
J = besselj(m, u);
dJ = 0.5*(besselj(m-1, u) - besselj(m+1, u));
JoverRho = zeros(size(rho));
awayFromAxis = rho ~= 0;
JoverRho(awayFromAxis) = J(awayFromAxis)./rho(awayFromAxis);
% m*J_m(q*rho)/rho has a finite limit for every integer m.
mJoverRho = m*JoverRho;
if abs(m) == 1
    mJoverRho(~awayFromAxis) = beam.q/2;
end

%% 3. Gas amplitudes -- retain the original vector-field equations.
Er = -(beam.omega*beam.muGas/beam.q^2)*beam.Ch*mJoverRho ...
    + 1i*(beam.betaGas/beam.q)*beam.Ce*dJ;
Ep = -(beam.betaGas/beam.q^2)*beam.Ce*mJoverRho ...
    - 1i*(beam.omega*beam.muGas/beam.q)*beam.Ch*dJ;
Ez = beam.Ce*J;
Ex = Er.*cos(phi) - Ep.*sin(phi);
Ey = Er.*sin(phi) + Ep.*cos(phi);

%% 4. Liquid projection -- preserve the source's returned lE components.
% The old 6x6 inverse and reflected-field calculations did not contribute
% to these liquid outputs or either requested figure, so they are omitted.
if geometryMode == "legacy"
    % These meshgrid axes reproduce X/Y as built late in the old function.
    [Xt, Yt] = meshgrid(grid.xM, grid.yM);
    % Compatibility: the old loop overwrote i with the number of x samples.
    ellipseRatio = sqrt(1 + (1-1/beam.indexRatio^2)*tan(numel(grid.xM))^2);
else
    Xt = X;
    Yt = Y;
    ellipseRatio = sqrt(1 + (1-1/beam.indexRatio^2)*tan(beam.incidenceRad)^2);
end
tangentAngle = atan(-Xt./(Yt*ellipseRatio^2));
tangentAngle(tangentAngle <= 0) = tangentAngle(tangentAngle <= 0) + pi;
tangentAngle(Xt == 0 & Yt == 0) = 0;
field.gas = struct('rho', Er, 'phi', Ep, 'z', Ez);
field.liquid = struct( ...
    'rho', cos(phi).*Ex + sin(phi).*Ey, ...
    'phi', cos(tangentAngle).*Ex + sin(tangentAngle).*Ey, ...
    'z', Ez);
field.phi = phi;
field.axialM = axialM;
end
