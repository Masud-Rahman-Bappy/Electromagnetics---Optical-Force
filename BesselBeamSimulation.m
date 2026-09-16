classdef BesselBeamSimulation < handle
    properties (SetAccess = private)
        Config (1,1) struct
        Results (1,1) struct = struct()
    end
    properties (Access = private)
        Grid (1,1) struct
        Beams (1,2) struct = repmat(struct(),1,2)
        Metrics (1,1) struct
        SpatialCache (1,2) cell = {[], []}
        PreparationSeconds (1,1) double = 0
    end

    methods
        function obj = BesselBeamSimulation(config)
            % Block 1: validate parameters and prepare immutable run constants.
            arguments
                config (1,1) struct
            end
            preparationTimer = tic;
            BesselBeamSimulation.validateConfig(config);
            obj.Config = config;
            obj.Grid = BesselBeamSimulation.makeGrid(config.grid);
            [obj.Beams, obj.Metrics] = BesselBeamSimulation.prepareBeams(config);
            % Spatial amplitudes are z-independent for the default equal angles.
            for b = 1:2
                if abs(obj.Beams(b).tiltRad) < 1e-14
                    obj.SpatialCache{b} = BesselBeamGasLiquidForceFunciton( ...
                        obj.Beams(b), obj.Grid, 0, config.numerics.geometryMode);
                end
            end
            obj.PreparationSeconds = toc(preparationTimer);
        end

        function result = run(obj)
            % Block 2: stream through z; only intensity volumes are retained.
            c = obj.Config;
            shape = [numel(obj.Grid.xM), numel(obj.Grid.yM), numel(obj.Grid.zM)];
            raw = zeros(shape, c.numerics.storageClass);
            normalized = zeros(shape, c.numerics.storageClass);
            gas = [];
            if c.numerics.storeGasIntensity
                gas = zeros(shape, c.numerics.storageClass);
            end
            timer = tic;
            evaluations = sum(~cellfun(@isempty, obj.SpatialCache));
            for iz = 1:shape(3)
                zM = obj.Grid.zM(iz);
                totalL = complex(zeros(shape(1), shape(2), 3));
                if c.numerics.storeGasIntensity
                    totalG = complex(zeros(shape(1), shape(2), 3));
                end
                for b = 1:2
                    beam = obj.Beams(b);
                    f = obj.SpatialCache{b};
                    if isempty(f)
                        f = BesselBeamGasLiquidForceFunciton( ...
                            beam, obj.Grid, zM, c.numerics.geometryMode);
                        evaluations = evaluations + 1;
                    else
                        f.axialM = zM + zeros(shape(1),shape(2));
                        if c.numerics.geometryMode == "legacy"
                            f.axialM(obj.Grid.X == 0 & obj.Grid.Y == 0) = 0;
                        end
                    end
                    % Remove the common optical carrier; only beat phase affects |E|^2.
                    phase = exp(1i*(beam.betaLiquid*f.axialM ...
                        - beam.deltaOmega*c.timeS + beam.order*f.phi));
                    totalL(:,:,1) = totalL(:,:,1) + f.liquid.rho.*phase;
                    totalL(:,:,2) = totalL(:,:,2) + f.liquid.phi.*phase;
                    totalL(:,:,3) = totalL(:,:,3) + f.liquid.z.*phase;
                    if c.numerics.storeGasIntensity
                        phaseG = exp(1i*(beam.betaGas*f.axialM ...
                            - beam.deltaOmega*c.timeS + beam.order*f.phi));
                        totalG(:,:,1) = totalG(:,:,1) + f.gas.rho.*phaseG;
                        totalG(:,:,2) = totalG(:,:,2) + f.gas.phi.*phaseG;
                        totalG(:,:,3) = totalG(:,:,3) + f.gas.z.*phaseG;
                    end
                end
                plane = sum(abs(totalL).^2, 3);
                peak = max(plane, [], 'all');
                if ~all(isfinite(plane), 'all') || peak <= 0
                    error('Bessel:InvalidIntensity', ...
                        'Nonfinite or zero intensity at axial plane %d.', iz);
                end
                raw(:,:,iz) = plane;
                normalized(:,:,iz) = plane/peak;
                if c.numerics.storeGasIntensity
                    gas(:,:,iz) = sum(abs(totalG).^2, 3);
                end
            end
            % Block 3: package data and metadata in one documented structure.
            result.config = c;
            result.grid = struct('xM', obj.Grid.xM*obj.Beams(1).stretch, ...
                'yM', obj.Grid.yM, 'zM', obj.Grid.zM);
            result.intensity.rawLiquid = raw; % Original electric-field sum |E|^2.
            result.intensity.normalizedLiquid = normalized; % Each z plane has max=1.
            result.intensity.rawGas = gas;
            result.metrics = obj.Metrics;
            result.performance.simulationSeconds = toc(timer);
            result.performance.preparationSeconds = obj.PreparationSeconds;
            result.performance.totalComputeSeconds = obj.PreparationSeconds ...
                + result.performance.simulationSeconds;
            result.performance.spatialEvaluations = evaluations;
            result.metadata.geometryMode = c.numerics.geometryMode;
            result.metadata.createdAt = datetime('now');
            result.metadata.arrayOrder = "x-index, y-index, z-index";
            result.metadata.model = "Original liquid-field projection";
            obj.Results = result;
            if c.numerics.verbose
                fprintf('%s: %d x %d x %d, %.3f s, %d spatial evaluations.\n', ...
                    c.name, shape(1), shape(2), shape(3), ...
                    result.performance.totalComputeSeconds, evaluations);
                fprintf('  Pattern vz = %.6g m/s; omega = %.6g rad/s.\n', ...
                    obj.Metrics.axialVelocityMps, obj.Metrics.angularVelocityRadps);
            end
        end

        function fig = plot(obj)
            % Block 4: reproduce each preset's slices, isosurface, and marker.
            if ~isfield(obj.Results, 'intensity')
                error('Bessel:RunFirst','Call run() before plot().');
            end
            r = obj.Results;
            c = obj.Config;
            % Use micrometers for readable axes; preserve original display-axis swap.
            x = r.grid.xM*1e6;
            y = r.grid.yM*1e6;
            z = r.grid.zM*1e6;
            V = double(r.intensity.normalizedLiquid);
            % slice interpolation is most robust with increasing axes.
            [y, permutation] = sort(y);
            V = V(:,permutation,:);
            [H, W, Z] = meshgrid(y, x, z);
            fig = figure('Name',char(c.name),'NumberTitle','off', ...
                'Color','white','Visible',char(c.plot.visible));
            ax = axes('Parent',fig);
            hold(ax,'on');
            if c.plot.kind == "adp"
                obj.drawSlices(ax,H,W,Z,V,y(end),x(end),[z(1),z(end),0]);
                [px,py] = BesselBeamSimulation.pickPeak(r);
                pz = 0;
                % A vertical plane through the selected transverse peak.
                a = atan2(py,px);
                s = linspace(-hypot(max(abs(x)),max(abs(y))), ...
                    hypot(max(abs(x)),max(abs(y))),c.plot.sliceResolution);
                [S,T] = meshgrid(s,linspace(z(1),z(end),c.plot.sliceResolution));
                obj.drawSlices(ax,H,W,Z,V,S*sin(a),S*cos(a),T);
                px = px*1e6;
                py = py*1e6;
            else
                % Preserve the source's illustrative projected-speed marker.
                pz = abs(r.metrics.legacyProjectedVelocityMps)*c.timeS*1e6;
                if ~isfinite(pz)
                    pz = 0;
                end
                requestedZ = [z(1), pz];
                requestedZ = requestedZ(requestedZ >= z(1) & requestedZ <= z(end));
                obj.drawSlices(ax,H,W,Z,V,y(1),[x(1),0],requestedZ);
                px = c.particle.markerXYM(1)*1e6;
                py = c.particle.markerXYM(2)*1e6;
            end
            geometry = isosurface(H,W,Z,V,c.plot.isoLevel);
            if ~isempty(geometry.vertices)
                surface = patch(ax,geometry,'FaceColor',[0.98 0.70 0.16], ...
                    'EdgeColor','none','FaceAlpha',c.plot.isoAlpha);
                isonormals(H,W,Z,V,surface);
            end
            [sx,sy,sz] = sphere(c.plot.sphereResolution);
            radius = c.particle.displayRadiusM*1e6;
            if c.plot.kind == "adp"
                horizontal = py;
                vertical = px;
            else
                % The original helical script places xc on its first display axis.
                horizontal = px;
                vertical = py;
            end
            surf(ax,horizontal+radius*sx,vertical+radius*sy,pz+radius*sz, ...
                'FaceColor',[0.2 0.9 0.3],'EdgeColor','none','FaceLighting','gouraud');
            colormap(ax,parula(256));
            caxis(ax,[0 1]);
            cb = colorbar(ax);
            cb.Label.String = 'Plane-normalized |E|^2';
            xlabel(ax,'y display coordinate (\mum)');
            ylabel(ax,'x display coordinate (\mum)');
            zlabel(ax,'z (\mum)');
            title(ax,sprintf('%s | t = %.3g microseconds',c.name,c.timeS*1e6), ...
                'Interpreter','none');
            set(ax,'Color','black','FontSize',11);
            axis(ax,[min(y),max(y),min(x),max(x),min(z),max(z)]);
            daspect(ax,[1 1 1]);
            view(ax,c.plot.viewDirection);
            box(ax,'on');
            camlight(ax,'headlight');
            lighting(ax,'gouraud');
            hold(ax,'off');
        end

        function saveOutputs(obj, fig)
            % Block 5: optional reproducible exports; no files are saved by default.
            arguments
                obj (1,1) BesselBeamSimulation
                fig (1,1) matlab.ui.Figure
            end
            o = obj.Config.output;
            if ~(o.saveFigures || o.saveResults)
                return
            end
            if ~isfolder(o.folder)
                mkdir(o.folder);
            end
            stem = fullfile(o.folder,char(obj.Config.id));
            if o.saveFigures
                savefig(fig,[stem '.fig']);
                exportgraphics(fig.CurrentAxes,[stem '.png'],'Resolution',300);
            end
            if o.saveResults
                result = obj.Results; %#ok<NASGU>
                save([stem '.mat'],'result','-v7.3');
            end
        end
    end

    methods (Static)
        function c = preset(name)
            % Block 6: nested structures group physical, grid, and display settings.
            arguments
                name (1,1) string {mustBeMember(name,["ADP15kHz","Helical5kHz"])}
            end
            c.id = name;
            c.name = "ADP 15 kHz";
            c.timeS = 0;
            c.beam = struct('wavelengthM',1e-6,'coneAngleDeg',40, ...
                'orders',[2 -2],'incidenceDeg',[45 45], ...
                'electricAmplitude',1.5e6,'magneticAmplitude',0, ...
                'frequencyDifferenceHz',15000, ...
                'betaDifferenceFraction',0.16090036886470404, ...
                'ringRadiiM',[5.403245142485066e-7,4.802389364443375e-7], ...
                'direction',1);
            c.medium = struct('muRelativeGas',1,'epsilonRelativeGas',1, ...
                'muRelativeLiquid',1,'epsilonRelativeLiquid',1.7689);
            c.grid = struct('xHalfWidthM',1e-6,'yHalfWidthM',1e-6, ...
                'xStepM',20e-9,'yStepM',20e-9, ...
                'zLimitsM',[-3e-6 3e-6],'zStepM',60e-9);
            c.particle = struct('displayRadiusM',150e-9, ...
                'markerXYM',[-600e-9 0]);
            c.plot = struct('kind',"adp",'isoLevel',0.99,'isoAlpha',0.75, ...
                'viewDirection',[-1 -1 1],'sphereResolution',50, ...
                'sliceResolution',160,'visible',"on");
            c.numerics = struct('geometryMode',"legacy",'storageClass',"double", ...
                'storeGasIntensity',false,'verbose',true);
            c.output = struct('saveFigures',false,'saveResults',false,'folder','results');
            if name == "Helical5kHz"
                c.name = "Helical pattern 5 kHz";
                c.timeS = 135e-6;
                c.beam.electricAmplitude = 2.5e6;
                c.beam.frequencyDifferenceHz = 5000;
                c.beam.betaDifferenceFraction = 0.45;
                c.grid.xHalfWidthM = 1.5e-6;
                c.grid.yHalfWidthM = 1.5e-6;
                c.grid.xStepM = 15e-9;
                c.grid.yStepM = 15e-9;
                c.particle.displayRadiusM = 180e-9;
                c.plot.kind = "helical";
                c.plot.isoLevel = 0.98;
                c.plot.viewDirection = [1 1 1];
                c.plot.sphereResolution = 100;
            end
        end
    end

    methods (Static, Access = private)
        function validateConfig(c)
            % Block 7: reject inconsistent inputs before allocating large arrays.
            validateattributes(c.timeS,{'double'},{'scalar','real','finite'});
            validateattributes(c.beam.wavelengthM,{'double'},{'scalar','positive','finite'});
            validateattributes(c.beam.coneAngleDeg,{'double'},{'scalar','>',0,'<',90});
            validateattributes(c.beam.orders,{'double'},{'size',[1 2],'integer','finite'});
            validateattributes(c.beam.incidenceDeg,{'double'},{'size',[1 2],'real','finite','>=',0,'<',90});
            validateattributes(c.beam.ringRadiiM,{'double'},{'size',[1 2],'positive','finite'});
            validateattributes(c.beam.frequencyDifferenceHz,{'double'},{'scalar','nonnegative','finite'});
            validateattributes(c.beam.betaDifferenceFraction,{'double'},{'scalar','>=',0,'<',1});
            validateattributes(c.beam.electricAmplitude,{'double'},{'scalar','finite'});
            validateattributes(c.beam.magneticAmplitude,{'double'},{'scalar','finite'});
            validateattributes(c.beam.direction,{'double'},{'scalar','integer','finite'});
            assert(any(c.beam.direction == [-1 1]),'Bessel:Direction','direction must be +1 or -1.');
            assert(any(c.numerics.geometryMode == ["legacy","consistent"]), ...
                'Bessel:GeometryMode','geometryMode must be legacy or consistent.');
            assert(any(c.numerics.storageClass == ["double","single"]), ...
                'Bessel:StorageClass','storageClass must be double or single.');
            values = [c.medium.muRelativeGas,c.medium.epsilonRelativeGas, ...
                c.medium.muRelativeLiquid,c.medium.epsilonRelativeLiquid];
            validateattributes(values,{'double'},{'positive','finite','real'});
            values = [c.grid.xHalfWidthM,c.grid.yHalfWidthM,c.grid.xStepM, ...
                c.grid.yStepM,c.grid.zStepM,c.particle.displayRadiusM];
            validateattributes(values,{'double'},{'positive','finite','real'});
            validateattributes(c.grid.zLimitsM,{'double'},{'size',[1 2],'finite','real'});
            assert(c.grid.zLimitsM(2)>c.grid.zLimitsM(1),'Bessel:Grid','z limits must increase.');
            validateattributes(c.plot.isoLevel,{'double'},{'scalar','>',0,'<',1});
            if c.numerics.geometryMode == "legacy"
                assert(c.grid.xHalfWidthM == c.grid.yHalfWidthM && ...
                    c.grid.xStepM == c.grid.yStepM,'Bessel:LegacyGrid', ...
                    'Legacy mode requires identical x/y grids; use consistent mode for rectangular grids.');
                assert(c.beam.incidenceDeg(1) == c.beam.incidenceDeg(2), ...
                    'Bessel:UnequalIncidence', ...
                    'This figure workflow requires equal incidence angles; unequal liquid grids need resampling.');
            else
                % Different refracted transverse grids require resampling, not direct addition.
                assert(c.beam.incidenceDeg(1) == c.beam.incidenceDeg(2), ...
                    'Bessel:UnequalIncidence', ...
                    'Consistent mode currently requires equal incidence angles for a shared liquid grid.');
            end
        end

        function grid = makeGrid(c)
            % Block 8: force odd transverse sample counts and an exact center.
            nx = round(c.xHalfWidthM/c.xStepM);
            ny = round(c.yHalfWidthM/c.yStepM);
            nz = round(diff(c.zLimitsM)/c.zStepM);
            assert(nx>=1 && ny>=1 && nz>=2,'Bessel:GridSize','Use at least three samples per axis.');
            grid.xM = linspace(-c.xHalfWidthM,c.xHalfWidthM,2*nx+1);
            grid.yM = linspace(c.yHalfWidthM,-c.yHalfWidthM,2*ny+1);
            grid.xM(nx+1) = 0;
            grid.yM(ny+1) = 0;
            grid.zM = linspace(c.zLimitsM(1),c.zLimitsM(2),nz+1);
            [grid.X,grid.Y] = ndgrid(grid.xM,grid.yM);
        end

        function [beams, metrics] = prepareBeams(c)
            % Block 9: calculate frequencies, wavevectors, refraction, and tilts once.
            mu0 = 4*pi*1e-7;
            epsilon0 = 8.854187817e-12;
            muG = c.medium.muRelativeGas*mu0;
            epsG = c.medium.epsilonRelativeGas*epsilon0;
            nRatio = sqrt(c.medium.muRelativeLiquid*c.medium.epsilonRelativeLiquid ...
                /(c.medium.muRelativeGas*c.medium.epsilonRelativeGas));
            speed = 1/sqrt(muG*epsG);
            incidence = c.beam.incidenceDeg*pi/180;
            assert(all(abs(sin(incidence)/nRatio)<=1),'Bessel:Refraction', ...
                'The present propagating-wave model does not support total internal reflection.');
            refracted = asin(sin(incidence)/nRatio);
            if c.numerics.geometryMode == "legacy"
                angleForCos = c.beam.incidenceDeg; % Original degree/radian mismatch.
            else
                angleForCos = incidence;
            end
            separation = refracted(1)-c.beam.direction*refracted(2);
            denominator = c.beam.ringRadiiM(2)/c.beam.ringRadiiM(1) ...
                *(cos(refracted(2))/cos(angleForCos(2))) ...
                /(cos(refracted(1))/cos(angleForCos(1))) + cos(separation);
            if c.numerics.geometryMode == "legacy"
                tilt1 = atan(sin(separation)/denominator);
            else
                tilt1 = atan2(sin(separation),denominator);
            end
            tilts = [tilt1, -(separation-tilt1)];
            k1 = 2*pi/c.beam.wavelengthM;
            omega1 = speed*k1;
            % Preserve Delw=delf/(3e8/lambda), including the original c approximation.
            relativeDifference = c.beam.frequencyDifferenceHz/(3e8/c.beam.wavelengthM);
            deltaOmega = omega1*relativeDifference;
            omega = [omega1,omega1-deltaOmega];
            beta1 = k1*cosd(c.beam.coneAngleDeg);
            beta = beta1*[1,1-c.beam.betaDifferenceFraction];
            k = omega/speed;
            assert(all(k.^2>beta.^2),'Bessel:Wavevector','Both transverse wavevectors must be real and nonzero.');
            base = struct('order',0,'omega',0,'deltaOmega',0,'betaGas',0, ...
                'betaLiquid',0,'q',0,'Ce',c.beam.electricAmplitude, ...
                'Ch',c.beam.magneticAmplitude,'muGas',muG,'indexRatio',nRatio, ...
                'incidenceRad',0,'tiltRad',0,'stretch',0);
            beams = repmat(base,1,2);
            for b = 1:2
                beams(b).order = c.beam.orders(b);
                beams(b).omega = omega(b);
                beams(b).deltaOmega = -(b-1)*deltaOmega;
                beams(b).betaGas = beta(b);
                beams(b).betaLiquid = nRatio*beta(b);
                beams(b).q = sqrt(k(b)^2-beta(b)^2);
                beams(b).incidenceRad = incidence(b);
                beams(b).tiltRad = tilts(b);
                stretchIndex = nRatio;
                if c.numerics.geometryMode == "legacy"
                    stretchIndex = 1.33;
                end
                beams(b).stretch = sqrt(1+(1-1/stretchIndex^2)*tan(incidence(b))^2) ...
                    /cos(tilts(b));
            end
            deltaBeta = nRatio*(beta(1)*cos(tilts(1))-beta(2)*cos(tilts(2)));
            if abs(deltaBeta) < 1e-12
                vz = NaN; % No finite axial pattern velocity in the degenerate case.
            else
                vz = -deltaOmega/deltaBeta;
            end
            deltaOrder = c.beam.orders(1)-c.beam.orders(2);
            if deltaOrder == 0
                angular = NaN;
            else
                angular = deltaOmega/deltaOrder;
            end
            projected = NaN;
            if isfinite(vz) && isfinite(angular) && angular ~= 0
                % Retained ONLY to reproduce the old illustrative sphere location.
                % atan(vz/angular) is not dimensionally a physical helix angle.
                projected = vz*cos(atan(vz/angular))^2;
            end
            metrics = struct('axialVelocityMps',vz,'angularVelocityRadps',angular, ...
                'legacyProjectedVelocityMps',projected, ...
                'deltaOmegaRadps',deltaOmega,'deltaBetaPerM',deltaBeta, ...
                'actualFrequencyDifferenceHz',deltaOmega/(2*pi), ...
                'tiltRad',tilts,'refractedAnglesRad',refracted);
        end

        function [x,y] = pickPeak(result)
            % Block 10: choose the fourth near-maximum when present, else the maximum.
            [~,iz] = min(abs(result.grid.zM));
            plane = result.intensity.normalizedLiquid(:,:,iz);
            candidates = find(plane > 0.999);
            if isempty(candidates)
                [~,idx] = max(plane(:));
            else
                idx = candidates(min(4,numel(candidates)));
            end
            [ix,iy] = ind2sub(size(plane),idx);
            x = result.grid.xM(ix);
            y = result.grid.yM(iy);
        end

        function drawSlices(ax,X,Y,Z,V,xs,ys,zs)
            h = slice(ax,X,Y,Z,V,xs,ys,zs);
            set(h,'FaceColor','interp','EdgeColor','none');
        end
    end
end
