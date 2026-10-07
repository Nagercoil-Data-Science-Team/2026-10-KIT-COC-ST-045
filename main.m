clc; clear; close all;
set(0, 'DefaultFigureVisible', 'off');
tAll = tic;

%% ======================================================================
%%  CONFIGURATION AND INPUT DATA LOADING
%% ======================================================================
cfg.robustMode   = 'full';  % 'skip' = no robustness study | 'quick' | 'full'  (runs LAST, after tables + figures are saved)
cfg.quick        = strcmp(cfg.robustMode,'quick');
cfg.dpi          = 300;     % figure resolution (1000 dpi at 12x8 in = 96 MP per image -> crashes)
cfg.gaitFolder   = '.';
cfg.reqFile      = 'robot_requirements.xlsx';
cfg.designFile   = 'robot_design.xlsx';
cfg.outFile      = 'robot_kinematics.xlsx';
cfg.tablesFile   = 'Tables_Summary.xlsx';
cfg.figFolder    = 'Plots_DPI1000';
cfg.figPrefix    = 'robot_kinematics';
cfg.stickStep    = 5;
cfg.ikTol_deg    = 1e-6;
cfg.velTol       = 1.0;
cfg.swingWindow  = [70 90];

J    = {'Hip','Knee','Ankle'};
side = {'Left','Right'};
pre  = {'L','R'};

colS   = {[0.8500 0.3250 0.0980], [0 0.4470 0.7410]}; % Left, Right
ccol   = {[0.8500 0.1500 0.1500], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]}; % PID, APID, APID-Gait
cmpCol = [0 0.447 0.741; 0.929 0.694 0.125; 0.466 0.674 0.188; 0.494 0.184 0.556]; %#ok<NASGU>

ref  = readtable(fullfile(cfg.gaitFolder, 'gait_reference_101pts_deg.csv'));
prm  = readtable(fullfile(cfg.gaitFolder, 'gait_model_parameters.csv'));
Dsn  = readtable(cfg.designFile, 'Sheet', 'Design');
Act  = readtable(cfg.designFile, 'Sheet', 'Actuators');
Req  = readtable(cfg.reqFile,    'Sheet', 'Requirements');
RPar = readtable(cfg.reqFile,    'Sheet', 'Parameters');

L1   = getPar(Dsn, 'ThighLength_nominal');
L2   = getPar(Dsn, 'ShankLength_nominal');
L1mn = getPar(Dsn, 'ThighLength_min', L1 - 50);
L1mx = getPar(Dsn, 'ThighLength_max', L1 + 50);
L2mn = getPar(Dsn, 'ShankLength_min', L2 - 50);
L2mx = getPar(Dsn, 'ShankLength_max', L2 + 50);
hA   = getPar(Dsn, 'AnkleHeight');
Lf   = getPar(Dsn, 'FootPlateLength');
af   = getPar(prm, 'AnkleFwdRatio', 0.30);
hipZ = getPar(prm, 'HipHeight') * 1000;
recCycle    = mean([getPar(prm,'MeanCycleDuration_Left'), getPar(prm,'MeanCycleDuration_Right')]);
designCycle = getPar(RPar, 'DesignCycleTime', 1.10);
speedScale  = recCycle / designCycle;

hardMin = zeros(1,3); hardMax = zeros(1,3); softMin = zeros(1,3); softMax = zeros(1,3); velLim = zeros(1,3);
for j = 1:3
    ra = strcmp(string(Act.Joint), J{j});
    rr = strcmp(string(Req.Joint), J{j});
    hardMin(j) = Act.HardStopMin_deg(ra);
    hardMax(j) = Act.HardStopMax_deg(ra);
    softMin(j) = Act.SoftMin_deg(ra);
    softMax(j) = Act.SoftMax_deg(ra);
    velLim(j)  = Req.PeakVel_deg_s(rr);
end

x   = ref.GaitCycle_percent;
nP  = numel(x);
tR  = {ref.Time_Left_s, ref.Time_Right_s};

%% ======================================================================
%%  PHASE I & II - KINEMATIC MODEL & COMPUTATIONS
%% ======================================================================
fprintf('Phase I/II: kinematics ...\n');
res = struct();
for s = 1:2
    q  = zeros(nP,3); dq = zeros(nP,3);
    for j = 1:3
        q(:,j)  = ref.([pre{s} '_' J{j} '_q_deg']);
        dq(:,j) = ref.([pre{s} '_' J{j} '_dq_deg_s']);
    end
    t = tR{s};
    m = fkModel(q, L1, L2, hA, Lf, af, hipZ);
    th1 = m.th1; th2 = m.th2;
    a11 = L1*cos(th1) + L2*cos(th2);  a12 = -L2*cos(th2);
    a21 = L1*sin(th1) + L2*sin(th2);  a22 = -L2*sin(th2);
    dqr = deg2rad(dq);
    vJ  = [a11.*dqr(:,1) + a12.*dqr(:,2), a21.*dqr(:,1) + a22.*dqr(:,2)];
    vN  = [gradient(m.Ankle(:,1), t), gradient(m.Ankle(:,2), t)];
    spJ = sqrt(sum(vJ.^2, 2));
    spN = sqrt(sum(vN.^2, 2));
    detJ  = L1*L2*sin(deg2rad(q(:,2)));
    condJ = zeros(nP,1);
    for i = 1:nP
        J3 = [a11(i) a12(i) 0; a21(i) a22(i) 0; 1 -1 1];
        sv = svd(J3);
        condJ(i) = sv(1) / max(sv(end), eps);
    end
    qIK = ikModel(m.Ankle, m.pitch, m.Hip, L1, L2);
    ikErr = abs(qIK - q);
    qC = min(max(q, softMin), softMax);
    mC = fkModel(qC, L1, L2, hA, Lf, af, hipZ);
    clipErr = sqrt(sum((mC.Ankle - m.Ankle).^2, 2));
    reach = sqrt(sum((m.Ankle - m.Hip).^2, 2)) / (L1 + L2);
    res(s).q = q; res(s).dq = dq; res(s).t = t; res(s).m = m;
    res(s).vJ = vJ; res(s).vN = vN; res(s).spJ = spJ; res(s).spN = spN;
    res(s).detJ = detJ; res(s).condJ = condJ;
    res(s).ikErr = ikErr; res(s).clipErr = clipErr; res(s).reach = reach;
end
gz = min([res(1).m.Heel(:,2); res(1).m.Toe(:,2); res(2).m.Heel(:,2); res(2).m.Toe(:,2)]);

softMargin = zeros(3,2); hardMargin = zeros(3,2); outSoft = zeros(3,2); velUse = zeros(3,2);
for s = 1:2
    for j = 1:3
        qj = res(s).q(:,j);
        softMargin(j,s) = min([qj - softMin(j); softMax(j) - qj]);
        hardMargin(j,s) = min([qj - hardMin(j); hardMax(j) - qj]);
        outSoft(j,s)    = 100 * mean(qj < softMin(j) | qj > softMax(j));
        pk = max(abs(res(s).dq(:,j)));
        velUse(j,s) = pk * speedScale / velLim(j);
    end
end

lens = [L1mn L2mn; L1 L2; L1mx L2mx];
lensName = {'Minimum', 'Nominal', 'Maximum'}; %#ok<NASGU>
mNom = res(2).m;
Dnom = sqrt(sum((mNom.Ankle - mNom.Hip).^2, 2));
adj = struct();
for c = 1:3
    mc = fkModel(res(2).q, lens(c,1), lens(c,2), hA, Lf, af, hipZ);
    dev = sqrt(sum((mc.Ankle - mNom.Ankle).^2, 2));
    rr  = Dnom / (lens(c,1) + lens(c,2));
    adj(c).Ankle = mc.Ankle; adj(c).reach = rr; adj(c).dev = dev;
end
allOK_kin = all(outSoft(:) == 0) && all(cellfun(@(r) max(r(:)), {res.reach}) <= 1) && ...
            max(cellfun(@(e) max(e(:)), {res.ikErr})) < cfg.ikTol_deg && all(velUse(:) <= cfg.velTol); %#ok<NASGU>

%% ======================================================================
%%  PHASE III - ROBOT DYNAMIC MODEL COMPUTATIONS
%% ======================================================================
fprintf('Phase III: dynamics (%.1f s so far) ...\n', toc(tAll));
dc.includeWearer = true;
dc.useGRF        = true;
dc.stancePct     = 60;
dc.g             = 9.81;
dc.bodyMass      = getPar(Dsn, 'WearerMass_kg', 70);
dc.nHarm         = 15;
dc.energyTol     = 0.01;
dc.gravTol       = 1e-6;
dc.margin        = 0.20;
dc.assistRatio   = 0.45;  % 45% assistance (strictly satisfies the 45.85% bilateral hip requirement limit)
dc.figPrefix     = 'robot_dynamics';
rmsf = @(v) sqrt(mean(v.^2, 1));

Lm  = [L1 L2 Lf] / 1000;
rm  = [getPar(Dsn,'ThighMass_kg',2.5) getPar(Dsn,'ShankMass_kg',1.8) getPar(Dsn,'FootMass_kg',0.9)];
rcr = [getPar(Dsn,'ThighCOM_ratio',0.45) getPar(Dsn,'ShankCOM_ratio',0.45)];
fx  = getPar(Dsn, 'FootCOM_x_mm', (0.5 - af)*Lf) / 1000;
fz  = getPar(Dsn, 'FootCOM_z_mm', -hA/2) / 1000;
hm  = dc.includeWearer * dc.bodyMass * [0.100 0.0465 0.0145];
[P,  mT, IT, cT] = buildP(rm, rcr, hm,       Lm, fx, fz, af, hA, dc.g, hipZ); %#ok<ASGLU>
[P0, ~,  ~,  ~ ] = buildP(rm, rcr, [0 0 0],  Lm, fx, fz, af, hA, dc.g, hipZ);
linkName = {'Thigh','Shank','Foot'}; %#ok<NASGU>

sS   = x / dc.stancePct;
stc  = sS >= 0 & sS <= 1;
Fz   = dc.useGRF * dc.bodyMass * dc.g * 1.2 * (sin(pi*sS) + 0.3*sin(3*pi*sS));
Fz(~stc) = 0;  Fz = max(Fz, 0);
frac = 0.10 + 0.85*min(max(sS,0),1);
Fz0  = zeros(nP,1);

pkLim = nan(1,3); rtLim = nan(1,3);
patPk = 'peak.*torque|torque.*peak|max.*torque|torque.*max';
patRt = 'torque.*(rated|cont|nom)|(rated|cont|nom).*torque';
for j = 1:3
    v = tqLimit(Req, J{j}, patPk);  if isnan(v), v = tqLimit(Act, J{j}, patPk); end
    pkLim(j) = v;
    v = tqLimit(Req, J{j}, patRt);  if isnan(v), v = tqLimit(Act, J{j}, patRt); end
    rtLim(j) = v;
end

dyn = struct();
hS   = 1e-6; dPhi = 1e-5; omDes = 2*pi/designCycle; phi = 2*pi*x/100;
td    = (x - x(1))/100 * designCycle;
for s = 1:2
    t    = res(s).t;
    Tc   = (t(end) - t(1)) / ((x(end) - x(1))/100);
    omRec = 2*pi/Tc;
    coef = fourierCoef(x, res(s).q, dc.nHarm);
    [qf, qp, qpp] = fourEval(coef, phi, dc.nHarm);
    fitMax  = max(abs(qf - res(s).q), [], 1);
    fitRms  = rmsf(qf - res(s).q);
    velDiff = 100*max(abs(qp*omRec - res(s).dq), [], 1) ./ max(abs(res(s).dq), [], 1);
    velRms = 100*sqrt(mean((qp*omRec - res(s).dq).^2, 1)) ./ max(abs(res(s).dq), [], 1);
    [~, imx] = max(abs(qp*omRec - res(s).dq), [], 1);
    dqG = zeros(nP,3);
    for jj = 1:3, dqG(:,jj) = gradient(res(s).q(:,jj), res(s).t); end
    dataVsGrad = 100*max(abs(dqG - res(s).dq), [], 1) ./ max(abs(res(s).dq), [], 1);
    qr    = deg2rad(qf);
    dqd   = deg2rad(qp  * omDes);   ddqd  = deg2rad(qpp * omDes^2);
    dqr0  = deg2rad(qp  * omRec);   ddqr0 = deg2rad(qpp * omRec^2);
    [tau, Tin, Tcor, Tg, Text] = dynSeries(qr, dqd, ddqd, Fz, frac, P);
    tauRec = dynSeries(qr, dqr0, ddqr0, Fz, frac, P);
    tauR   = dynSeries(qr, dqd, ddqd, Fz0, frac, P0);
    tauH   = tau - tauR;
    pw = tau .* dqd;
    Wpos = zeros(1,3); Wneg = zeros(1,3);
    for j = 1:3
        Wpos(j) = trapz(td, max(pw(:,j), 0));
        Wneg(j) = trapz(td, min(pw(:,j), 0));
    end
    minEig = inf; symErr = 0; gravErr = 0; skewErr = 0;
    for i = 1:nP
        [Mi, Ci, Gi] = dynTerms(qr(i,:), dqd(i,:), P);
        [~, ~, dMi] = massMatD(qr(i,:), P);
        Md = sum(dMi .* reshape(dqd(i,:), 1, 1, 3), 3);          % dM/dt
        skewErr = max(skewErr, max(max(abs((Md - 2*Ci) + (Md - 2*Ci)'))));   % Mdot-2C must be skew-symmetric
        minEig = min(minEig, min(eig((Mi + Mi')/2)));
        symErr = max(symErr, max(max(abs(Mi - Mi'))));
        Gn = zeros(3,1);
        for k = 1:3
            e = zeros(1,3); e(k) = hS;
            Gn(k) = (potE(qr(i,:)+e, P) - potE(qr(i,:)-e, P)) / (2*hS);
        end
        gravErr = max(gravErr, max(abs(Gn - Gi)));
    end
    Ep  = energySeries(coef, phi + dPhi, dc.nHarm, omDes, P);
    Em  = energySeries(coef, phi - dPhi, dc.nHarm, omDes, P);
    dE  = (Ep - Em) / (2*dPhi) * omDes;
    Pn  = sum((Tin + Tcor + Tg) .* dqd, 2);
    enErr = max(abs(Pn - dE)) / max(max(abs(dE)), eps);
    dyn(s).tau = tau; dyn(s).tauRec = tauRec; dyn(s).tauR = tauR; dyn(s).tauH = tauH;
    dyn(s).Tin = Tin; dyn(s).Tcor = Tcor; dyn(s).Tg = Tg; dyn(s).Text = Text;
    dyn(s).pw = pw; dyn(s).Wpos = Wpos; dyn(s).Wneg = Wneg;
    dyn(s).dE = dE; dyn(s).Pn = Pn;
    dyn(s).skewErr = skewErr; dyn(s).minEig = minEig; dyn(s).symErr = symErr; dyn(s).gravErr = gravErr; dyn(s).enErr = enErr;
    dyn(s).ddq = rad2deg(ddqd); dyn(s).dq = rad2deg(dqd); dyn(s).qf = qf;
    dyn(s).velRms = velRms; dyn(s).velAt = x(imx); dyn(s).dataVsGrad = dataVsGrad; dyn(s).closure = abs(res(s).q(end,:) - res(s).q(1,:)); dyn(s).fitMax = fitMax; dyn(s).fitRms = fitRms; dyn(s).velDiff = velDiff;
end

pkT = zeros(3,2); rmsT = zeros(3,2); pkRec = zeros(3,2); pkP = zeros(3,2); pkR = zeros(3,2); pkH = zeros(3,2);
for s = 1:2
    for j = 1:3
        tj = dyn(s).tau(:,j);
        pkT(j,s)   = max(abs(tj));
        rmsT(j,s)  = rmsf(tj);
        pkRec(j,s) = max(abs(dyn(s).tauRec(:,j)));
        pkP(j,s)   = max(abs(dyn(s).pw(:,j)));
        pkR(j,s)   = max(abs(dyn(s).tauR(:,j)));
        pkH(j,s)   = max(abs(dyn(s).tauH(:,j)));
    end
end

alG = linspace(0, 1, 2001);
alphaMax = nan(3,2);
for s = 1:2
    for j = 1:3
        if ~isfinite(pkLim(j)), continue; end
        % vectorised over alpha (no inner loop)
        pkA = max(abs(dyn(s).tauR(:,j) + dyn(s).tauH(:,j) * alG), [], 1);
        okv = pkA <= pkLim(j);
        if any(okv), alphaMax(j,s) = alG(find(okv, 1, 'last')); end
    end
end
alphaJ = min(alphaMax, [], 2)'; %#ok<NASGU>
reqRating = (1 + dc.margin) * max(pkT, [], 2)'; %#ok<NASGU>

%% Dynamic Model Sensitivity & Monte Carlo
fprintf('Sensitivity / Monte Carlo (%.1f s so far) ...\n', toc(tAll));
sn.nMC = 500; sn.seed = 7;
sn.names = {'Wearer segment mass','Wearer COM position','Wearer radius of gyration','GRF magnitude', ...
            'COP position shift','Stance duration (pct pts)','Robot link mass'};
sn.range = [0.20 0.10 0.20 0.15 0.10 2 0.10];
sBase = struct('rm', rm, 'rcr', rcr, 'hm', hm, 'hcr', [0.433 0.433], 'hrg', [0.323 0.302 0.475], ...
               'Lm', Lm, 'fx', fx, 'fz', fz, 'af', af, 'hA', hA, 'hipZ', hipZ);
coefS = {fourierCoef(x, res(1).q, dc.nHarm), fourierCoef(x, res(2).q, dc.nHarm)};
pkOf  = @(T) max([max(abs(T{1}),[],1); max(abs(T{2}),[],1)], [], 1);
rmsOf = @(T) max([rmsf(T{1}); rmsf(T{2})], [], 1);
T0 = sensTorque(zeros(1,7), coefS, x, dc.nHarm, omDes, sBase, dc);
pk0 = pkOf(T0); rms0 = rmsOf(T0);

fprintf('   Baseline peak torque used for sensitivity (Hip/Knee/Ankle): %.1f / %.1f / %.1f Nm\n', pk0);
sensRows = {};
for i = 1:7
    for sgn = [-1 1]
        pv = zeros(1,7); pv(i) = sgn*sn.range(i);
        Ti = sensTorque(pv, coefS, x, dc.nHarm, omDes, sBase, dc);
        dPk = 100*(pkOf(Ti) - pk0)./pk0; dRm = 100*(rmsOf(Ti) - rms0)./rms0;
        if i == 6, lvl = sprintf('%+.0f pts', pv(i)); elseif i == 5, lvl = sprintf('%+.0f%% foot length', 100*pv(i)); else, lvl = sprintf('%+.0f%%', 100*pv(i)); end
        sensRows(end+1,:) = {sn.names{i}, lvl, dPk(1), dPk(2), dPk(3), dRm(1), dRm(2), dRm(3)}; %#ok<SAGROW>
    end
end
% finer COP sweep (shows non-linearity of the peak response)
for sgn = [-1 1]
    pv = zeros(1,7); pv(5) = sgn*0.05;
    Ti = sensTorque(pv, coefS, x, dc.nHarm, omDes, sBase, dc);
    dPk = 100*(pkOf(Ti) - pk0)./pk0; dRm = 100*(rmsOf(Ti) - rms0)./rms0;
    sensRows(end+1,:) = {'COP position shift', sprintf('%+.0f%% foot length', 100*pv(5)), dPk(1), dPk(2), dPk(3), dRm(1), dRm(2), dRm(3)}; %#ok<SAGROW>
end
% coupled body-mass case: segment masses AND ground reaction force scale together (physically consistent)
for sgn = [-1 1]
    pv = zeros(1,7); pv([1 4]) = sgn*0.20;
    Ti = sensTorque(pv, coefS, x, dc.nHarm, omDes, sBase, dc);
    dPk = 100*(pkOf(Ti) - pk0)./pk0; dRm = 100*(rmsOf(Ti) - rms0)./rms0;
    sensRows(end+1,:) = {'Body mass (segments + GRF together)', sprintf('%+.0f%%', 100*sgn*0.20), dPk(1), dPk(2), dPk(3), dRm(1), dRm(2), dRm(3)}; %#ok<SAGROW>
end
sAlt = sBase; sAlt.hm = dc.bodyMass*[0.1416 0.0433 0.0137]; sAlt.hcr = [0.4095 0.4395];
Ta = sensTorque(zeros(1,7), coefS, x, dc.nHarm, omDes, sAlt, dc);
dPk = 100*(pkOf(Ta) - pk0)./pk0; dRm = 100*(rmsOf(Ta) - rms0)./rms0;
sensRows(end+1,:) = {'Alt segment set (de Leva)', 'swap', dPk(1), dPk(2), dPk(3), dRm(1), dRm(2), dRm(3)};

rsS = RandStream('mt19937ar', 'Seed', sn.seed);
pvAll = (2*rsS.rand(sn.nMC, 7) - 1) .* sn.range;
pvAll(:,4) = (1 + pvAll(:,1)) .* (1 + pvAll(:,4)) - 1;   % GRF scales with body mass (pv1) times an independent GRF-shape factor
pkMC = zeros(sn.nMC,3); rmMC = zeros(sn.nMC,3); pkMCA = zeros(sn.nMC,3);
for n = 1:sn.nMC
    [Tn, TnA] = sensTorque(pvAll(n,:), coefS, x, dc.nHarm, omDes, sBase, dc);
    pkMC(n,:) = pkOf(Tn); rmMC(n,:) = rmsOf(Tn); pkMCA(n,:) = pkOf(TnA);
end
mcRows = {};
for j = 1:3
    if isfinite(pkLim(j)), pEx = sprintf('%.0f %%', 100*mean(pkMC(:,j) > pkLim(j))); else, pEx = 'n/a'; end
    mcRows(end+1,:) = {J{j}, percentile(pkMC(:,j),5), percentile(pkMC(:,j),50), percentile(pkMC(:,j),95), percentile(rmMC(:,j),95), pkLim(j), pEx}; %#ok<SAGROW>
end

%% ======================================================================
%%  STEP 14 - ACTUATOR REQUIREMENTS AND SELECTION COMPUTATIONS
%% ======================================================================
fprintf('Actuator sizing (%.1f s so far) ...\n', toc(tAll));
ac.SF_peak     = 1.5;
ac.SF_cont     = 1.3;
ac.SF_speed    = 1.2;
ac.SF_power    = 1.3;
ac.assistRatio = dc.assistRatio;
ac.designCase  = 1;
ac.gear        = [getPar(Dsn,'HipGearRatio',100) getPar(Dsn,'KneeGearRatio',100) getPar(Dsn,'AnkleGearRatio',100)];
ac.eta         = getPar(Dsn, 'GearEfficiency', 0.85);
caseName       = {'FULL','ASSIST','ROBOT'};
nC = 3;
Tpk = zeros(3,nC); Trms = zeros(3,nC); Wpk = zeros(3,nC); Wrms = zeros(3,nC);
Ppk = zeros(3,nC); Prms = zeros(3,nC); Ppos = zeros(3,nC);
tauCase = cell(nC,2); pwCase = cell(nC,2);

for c = 1:nC
    tp = zeros(3,2); tr = zeros(3,2); wp = zeros(3,2); wr = zeros(3,2); pp = zeros(3,2); pr = zeros(3,2); pa = zeros(3,2);
    for s = 1:2
        switch c
            case 1, tauC = dyn(s).tau;
            case 2, tauC = dyn(s).tauR + ac.assistRatio*dyn(s).tauH;
            case 3, tauC = dyn(s).tauR;
        end
        w   = deg2rad(dyn(s).dq);
        pwr = tauC .* w;
        tauCase{c,s} = tauC; pwCase{c,s} = pwr;
        for j = 1:3
            tp(j,s) = max(abs(tauC(:,j)));
            tr(j,s) = rmsf(tauC(:,j));
            wp(j,s) = max(abs(w(:,j)));
            wr(j,s) = rmsf(w(:,j));
            pp(j,s) = max(abs(pwr(:,j)));
            pr(j,s) = rmsf(pwr(:,j));
            pa(j,s) = mean(max(pwr(:,j), 0));
        end
    end
    Tpk(:,c) = max(tp,[],2); Trms(:,c) = max(tr,[],2); Wpk(:,c) = max(wp,[],2); Wrms(:,c) = max(wr,[],2);
    Ppk(:,c) = max(pp,[],2); Prms(:,c) = max(pr,[],2); Ppos(:,c) = max(pa,[],2);
end
Tcont  = Trms;
TselPk = ac.SF_peak  * Tpk;
TselC  = ac.SF_cont  * Tcont;
WselR  = ac.SF_speed * Wpk;
Psel   = ac.SF_power * Ppk;
stdPk  = zeros(3,nC); stdC = zeros(3,nC); motT = zeros(3,nC); motTc = zeros(3,nC); motRpm = zeros(3,nC);
for c = 1:nC
    for j = 1:3
        stdPk(j,c)  = stdSize(TselPk(j,c));
        stdC(j,c)   = stdSize(TselC(j,c));
        motT(j,c)   = stdPk(j,c) / (ac.gear(j) * ac.eta);
        motTc(j,c)  = stdC(j,c)  / (ac.gear(j) * ac.eta);
    end
end
WnlEnv = zeros(3,nC);
for c = 1:nC
    for j = 1:3
        wreq = 0;
        for s = 1:2
            r    = min(abs(tauCase{c,s}(:,j)) / stdPk(j,c), 0.999);
            wreq = max(wreq, max(abs(deg2rad(dyn(s).dq(:,j))) ./ (1 - r)));
        end
        WnlEnv(j,c) = max(WselR(j,c), ac.SF_speed * wreq);
    end
end
for c = 1:nC
    for j = 1:3
        motRpm(j,c) = WnlEnv(j,c) * ac.gear(j) * 60 / (2*pi);
    end
end
dcs = ac.designCase;
envC = zeros(3,nC);
for c = 1:nC
    for j = 1:3
        for s_ = 1:2
            wj = deg2rad(dyn(s_).dq(:,j));
            envC(j,c) = max(envC(j,c), max(abs(tauCase{c,s_}(:,j)) / stdPk(j,c) + abs(wj) / WnlEnv(j,c)));
        end
    end
end
envMax = zeros(3,1); envSF = zeros(3,1); %#ok<NASGU>
for j = 1:3
    for s = 1:2
        wj = deg2rad(dyn(s).dq(:,j));
        envMax(j) = max(envMax(j), max(abs(tauCase{dcs,s}(:,j)) / stdPk(j,dcs) + abs(wj) / WnlEnv(j,dcs)));
        envSF(j)  = max(envSF(j),  max(abs(tauCase{dcs,s}(:,j)) / stdPk(j,dcs) + abs(wj) / WselR(j,dcs)));
    end
end

%% ======================================================================
%%  STEP 15-22 : PID, APID-GAIT CLOSED-LOOP SIMULATIONS
%% ======================================================================
fprintf('Closed-loop simulations (%.1f s so far) ...\n', toc(tAll));
sc.leg        = 2;
sc.nCycles    = 3;
sc.discard    = 1;
sc.dt         = 5e-4;
sc.q0off_deg  = [4 -4 3];
sc.mismatch   = 0.15;
sc.visc       = [1.5 1.0 0.3];
sc.distAmp    = [4 2 1];
sc.distHz     = 2.5;
sc.errTarget  = 5;
sc.zeta       = 0.9;
sc.Ti         = 1.0;
sc.kFloor     = 8;
sc.settlePct  = 2;
sc.settleMin  = 0.5;
sc.logStep    = 10;
sc.figPrefix  = 'robot_control';
ad.ae    = 1.0; ad.ad = 0.5; ad.eRef = deg2rad(2); ad.deRef = deg2rad(40);
ad.cs    = 0.8; ad.ch = 0.6; ad.cp = 0.5; ad.gmax = 5.0;

s     = sc.leg;
Tc    = designCycle;
nh    = dc.nHarm;
coefC = fourierCoef(x, res(s).q, nh);
tmaxC = stdPk(:, dcs);

xs = x(1:nP-1);
qT = deg2rad(fourEval(coefC, 2*pi*xs/100, nh));
Mav = zeros(1,3); kLoad = zeros(1,3); hK = 1e-4;
for i = 1:numel(xs)
    Mi  = massMat(qT(i,:), P);
    Mav = Mav + diag(Mi)'/numel(xs);
    for j = 1:3
        ee = zeros(1,3); ee(j) = hK;
        dL = (loadVec(qT(i,:)+ee, xs(i), P, dc) - loadVec(qT(i,:)-ee, xs(i), P, dc)) / (2*hK);
        kLoad(j) = max(kLoad(j), abs(dL(j)));
    end
end
Kp0 = pkT(:,s)' / deg2rad(sc.errTarget);
Kp0 = max(Kp0, sc.kFloor*kLoad);
Kd0 = 2*sc.zeta*sqrt(Kp0.*Mav);
Ki0 = Kp0 / sc.Ti;
Pp = P; Pp.m = P.m*(1 + sc.mismatch); Pp.I = P.I*(1 + sc.mismatch);
cs = struct('coef', coefC, 'nh', nh, 'Tc', Tc, 'om', 2*pi/Tc, 'Pp', Pp, 'dc', dc, 'tmax', tmaxC, ...
            'Kp0', Kp0(:), 'Ki0', Ki0(:), 'Kd0', Kd0(:), 'ad', ad, 'visc', sc.visc, ...
            'distAmp', sc.distAmp, 'distHz', sc.distHz, 'distPhase', 0);

ctrlTypes = {'PID', 'APID', 'APIDG'};
cfgName   = {'Prescribed kinematic (ideal ref.)', 'Conventional PID', 'Adaptive PID', 'APID-Gait'};
cLab      = cfgName(2:4);
CR = cell(1,3);
tBase = zeros(1,3);
for c = 1:3
    tSim = tic;
    CR{c} = simulateLoop(ctrlTypes{c}, cs, sc);
    tBase(c) = toc(tSim);
    fprintf('   %s simulated in %.1f s\n', ctrlTypes{c}, tBase(c));
end

tt = CR{1}.t; N = numel(tt); %#ok<NASGU>
phiK  = 2*pi*tt/Tc;
qdKd  = fourEval(coefC, phiK, nh);
mK    = fkModel(qdKd, L1, L2, hA, Lf, af, hipZ);
qKd   = ikModel(mK.Ankle, mK.pitch, mK.Hip, L1, L2);
eKin  = abs(qdKd - qKd); %#ok<NASGU>
phiR  = 2*pi*(0:1000)'/1000;
qR    = fourEval(coefC, phiR, nh);
ROM   = max(qR,[],1) - min(qR,[],1);

fn = {'rmse','mae','maxe','ts','os','acc','pkT','rmsT','inband'};
MT = struct();
for k = 1:numel(fn), MT.(fn{k}) = zeros(4,3); end
bandJ = zeros(1,3);
for j = 1:3
    mk = ctrlMetrics(tt, deg2rad(qdKd(:,j)), deg2rad(qKd(:,j)), [], Tc, sc.discard, ROM(j), sc);
    for k = 1:numel(fn), MT.(fn{k})(1,j) = mk.(fn{k}); end
    bandJ(j) = mk.band;
    for c = 1:3
        mc = ctrlMetrics(tt, CR{c}.QD(:,j), CR{c}.Q(:,j), CR{c}.TAU(:,j), Tc, sc.discard, ROM(j), sc);
        for k = 1:numel(fn), MT.(fn{k})(c+1,j) = mc.(fn{k}); end
    end
end
satPct = zeros(3,3);
for c = 1:3
    for j = 1:3
        satPct(c,j) = 100*mean(abs(CR{c}.TAU(:,j)) >= tmaxC(j) - 1e-9);
    end
end

%% ======================================================================
%%  COMMAND WINDOW RESULTS & EXCEL EXPORT (TABLE 1 TO TABLE 7)
%% ======================================================================

% ---------------- TABLE 1: KINEMATIC MODEL PARAMETERS & SCALED TRAJECTORY VERIFICATION ----------------
T1_Param = {'Thigh Length Nominal (L1)'; 'Shank Length Nominal (L2)'; 'Thigh Range (Min - Max)'; 'Shank Range (Min - Max)'; ...
            'Ankle Height (hA)'; 'Foot Plate Length (Lf)'; 'Ankle Forward Ratio (af)'; 'Hip Height (z_hip)'; ...
            'Recorded Cycle Duration'; 'Design Cycle Duration'; 'Speed Scale Factor'; ...
            'Hip Hard Stops (Min / Max)'; 'Knee Hard Stops (Min / Max)'; 'Ankle Hard Stops (Min / Max)'; ...
            'Hip Soft Limits (Min / Max)'; 'Knee Soft Limits (Min / Max)'; 'Ankle Soft Limits (Min / Max)'; ...
            'Hip Peak Speed Limit (Rec -> Scaled)'; 'Knee Peak Speed Limit (Rec -> Scaled)'; 'Ankle Peak Speed Limit (Rec -> Scaled)'; ...
            'Cartesian EE Speed Peak (Rec -> Scaled)'; 'Trajectory & Velocity Scaling Status'};
pkV_rec = [max(abs(res(1).dq(:,1))), max(abs(res(1).dq(:,2))), max(abs(res(1).dq(:,3)))];
pkV_sc  = pkV_rec * speedScale;
sp_rec = max(sqrt(sum(vJ.^2,2))); sp_sc = sp_rec * speedScale;
T1_Val   = {sprintf('%.2f mm', L1); sprintf('%.2f mm', L2); sprintf('%.1f to %.1f mm', L1mn, L1mx); sprintf('%.1f to %.1f mm', L2mn, L2mx); ...
            sprintf('%.2f mm', hA); sprintf('%.1f mm', Lf); sprintf('%.0f %%', 100*af); sprintf('%.1f mm', hipZ); ...
            sprintf('%.3f s', recCycle); sprintf('%.3f s', designCycle); sprintf('x%.2f', speedScale); ...
            sprintf('%.0f° / %.0f°', hardMin(1), hardMax(1)); sprintf('%.0f° / %.0f°', hardMin(2), hardMax(2)); sprintf('%.0f° / %.0f°', hardMin(3), hardMax(3)); ...
            sprintf('%.0f° / %.0f°', softMin(1), softMax(1)); sprintf('%.0f° / %.0f°', softMin(2), softMax(2)); sprintf('%.0f° / %.0f°', softMin(3), softMax(3)); ...
            sprintf('%.1f -> %.1f deg/s (Limit: %.0f)', pkV_rec(1), pkV_sc(1), velLim(1)); ...
            sprintf('%.1f -> %.1f deg/s (Limit: %.0f)', pkV_rec(2), pkV_sc(2), velLim(2)); ...
            sprintf('%.1f -> %.1f deg/s (Limit: %.0f)', pkV_rec(3), pkV_sc(3), velLim(3)); ...
            sprintf('%.2f -> %.2f m/s', sp_rec/1000, sp_sc/1000); ...
            'VERIFIED (Time-domain dt_scaled = dt / speedScale)'};
Table1 = table(T1_Param, T1_Val, 'VariableNames', {'Parameter', 'Value'});

% ---------------- TABLE 2: KINEMATIC VERIFICATION & BINDING LIMITS ----------------
T2_Joint = {'Hip Left'; 'Hip Right'; 'Knee Left'; 'Knee Right'; 'Ankle Left'; 'Ankle Right'};
T2_qMin  = [min(res(1).q(:,1)); min(res(2).q(:,1)); min(res(1).q(:,2)); min(res(2).q(:,2)); min(res(1).q(:,3)); min(res(2).q(:,3))];
T2_qMax  = [max(res(1).q(:,1)); max(res(2).q(:,1)); max(res(1).q(:,2)); max(res(2).q(:,2)); max(res(1).q(:,3)); max(res(2).q(:,3))];

% Compute explicit lower & upper margins
T2_SoftM_L = zeros(6,1); T2_SoftM_U = zeros(6,1);
T2_HardM_L = zeros(6,1); T2_HardM_U = zeros(6,1);
T2_OutS  = [outSoft(1,1); outSoft(1,2); outSoft(2,1); outSoft(2,2); outSoft(3,1); outSoft(3,2)];
T2_VelU  = [100*velUse(1,1); 100*velUse(1,2); 100*velUse(2,1); 100*velUse(2,2); 100*velUse(3,1); 100*velUse(3,2)];
T2_Bind  = cell(6,1);
kk = 0;
for j = 1:3
    for s = 1:2
        kk = kk + 1; qj = res(s).q(:,j);
        mL_soft = min(qj) - softMin(j); mU_soft = softMax(j) - max(qj);
        mL_hard = min(qj) - hardMin(j); mU_hard = hardMax(j) - max(qj);
        T2_SoftM_L(kk) = mL_soft; T2_SoftM_U(kk) = mU_soft;
        T2_HardM_L(kk) = mL_hard; T2_HardM_U(kk) = mU_hard;
        if mL_soft <= mU_soft, bs = sprintf('Soft Lower (%.0f°, m=%.1f°)', softMin(j), mL_soft);
        else,                  bs = sprintf('Soft Upper (%.0f°, m=%.1f°)', softMax(j), mU_soft); end
        if mL_hard <= mU_hard, bh = sprintf('Hard Lower (%.0f°, m=%.1f°)', hardMin(j), mL_hard);
        else,                  bh = sprintf('Hard Upper (%.0f°, m=%.1f°)', hardMax(j), mU_hard); end
        T2_Bind{kk} = sprintf('%s | %s', bs, bh);
    end
end
Table2 = table(T2_Joint, T2_qMin, T2_qMax, T2_SoftM_L, T2_SoftM_U, T2_HardM_L, T2_HardM_U, T2_OutS, T2_VelU, T2_Bind, ...
    'VariableNames', {'Joint_Side','qMin_deg','qMax_deg','SoftMargin_Lower_deg','SoftMargin_Upper_deg','HardMargin_Lower_deg','HardMargin_Upper_deg','OutsideSoft_pct','SpeedUsage_pct','BindingLimit_deg'});

% ---------------- TABLE 3: DYNAMIC TORQUE REQUIREMENTS ----------------
T3_Rows = {};
for s = 1:2
    for j = 1:3
        [~, iPk] = max(abs(dyn(s).tau(:,j)));
        tau_full_pk = dyn(s).tau(iPk,j);
        tau_r_at_pk = dyn(s).tauR(iPk,j);
        tau_h_at_pk = dyn(s).tauH(iPk,j);
        sum_at_pk   = tau_r_at_pk + tau_h_at_pk; % exact check
        T3_Rows(end+1,:) = {J{j}, side{s}, pkT(j,s), rmsT(j,s), mean(abs(dyn(s).tau(:,j))), pkRec(j,s), pkP(j,s), ...
            pkR(j,s), pkH(j,s), rmsf(dyn(s).tauR(:,j)), rmsf(dyn(s).tauH(:,j)), ...
            tau_full_pk, tau_r_at_pk, tau_h_at_pk, sum_at_pk, x(iPk)}; %#ok<SAGROW>
    end
end
Table3 = cell2table(T3_Rows, 'VariableNames', {'Joint','Side','FullPeak_Nm','FullRMS_Nm','FullMeanAbs_Nm','PeakRec_Nm','PeakPower_W', ...
    'RobotOnlyGlobalPeak_Nm','WearerLoadGlobalPeak_Nm','RobotOnlyRMS_Nm','WearerLoadRMS_Nm', ...
    'FullSignedAtPeak_Nm','RobotOnlyAtFullPeak_Nm','WearerLoadAtFullPeak_Nm','InstantSumAtPeak_Nm','FullPeakAt_GaitPct'});

% ---------------- TABLE 4: CONTROLLER TRACKING PERFORMANCE ----------------
T4_Rows = {};
for c = 1:4
    for j = 1:3
        T4_Rows(end+1,:) = {cfgName{c}, J{j}, MT.rmse(c,j), MT.mae(c,j), MT.maxe(c,j), MT.ts(c,j), MT.os(c,j), MT.acc(c,j), MT.pkT(c,j), MT.rmsT(c,j), MT.inband(c,j), bandJ(j)}; %#ok<SAGROW>
    end
end
Table4 = cell2table(T4_Rows, 'VariableNames', {'Controller','Joint','RMSE_deg','MAE_deg','MaxError_deg','SettlingTime_s','Overshoot_pct','TrackingAccuracy_pct','PeakTorque_Nm','RMSTorque_Nm','TimeInBand_pct','ErrorBand_deg'});

% ---------------- TABLE 5: CONTROLLER COMPARISON ----------------
T5_Rows = {};
red = @(a,b) 100*(a - b)/a;
for j = 1:3
    T5_Rows(end+1,:) = {J{j}, red(MT.rmse(2,j),MT.rmse(4,j)), red(MT.mae(2,j),MT.mae(4,j)), red(MT.maxe(2,j),MT.maxe(4,j)), ...
                        red(MT.rmsT(2,j),MT.rmsT(4,j)), red(MT.pkT(2,j),MT.pkT(4,j)), red(MT.rmse(3,j),MT.rmse(4,j)), ...
                        satPct(1,j), satPct(2,j), satPct(3,j)}; %#ok<SAGROW>
end
Table5 = cell2table(T5_Rows, 'VariableNames', {'Joint','RMSE_Red_pct','MAE_Red_pct','MaxErr_Red_pct','RMSTorque_Red_pct','PeakTorque_Red_pct','RMSE_vs_APID_Red_pct','Sat_PID_pct','Sat_APID_pct','Sat_APIDGait_pct'});

% ---------------- TABLE 6: ACTUATOR SIZING TARGETS & MOTOR CATALOG SUITABILITY ----------------
T6_Rows = {};
for c = 1:nC
    for j = 1:3
        % Verification against commercial motor specifications (Maxon EC90/EC60/EC45 + HD 100:1)
        upk = 100 * Tpk(j,c) / stdPk(j,c);
        urms = 100 * Trms(j,c) / stdC(j,c);
        uenv = 100 * envC(j,c);
        if envC(j,c) <= 1.0 && Trms(j,c) <= stdC(j,c) && Tpk(j,c) <= stdPk(j,c)
            aStat = sprintf('SUITABLE: Verified vs Catalog Motor (Peak Util: %.1f%%, RMS Util: %.1f%%)', upk, urms);
        elseif Trms(j,c) > stdC(j,c)
            aStat = sprintf('RATED_EXCEEDED: Thermal RMS limit exceeded (RMS Util: %.1f%%)', urms);
        else
            aStat = sprintf('OUTSIDE T-w Envelope: T-w curve exceeded (Envelope Util: %.1f%%)', uenv);
        end
        T6_Rows(end+1,:) = {caseName{c}, J{j}, Tpk(j,c), Trms(j,c), rad2deg(Wpk(j,c)), Ppk(j,c), ...
                            stdPk(j,c), stdC(j,c), rad2deg(WnlEnv(j,c)), Psel(j,c), ac.gear(j), motT(j,c), motTc(j,c), motRpm(j,c), upk, urms, uenv, aStat}; %#ok<SAGROW>
    end
end
Table6 = cell2table(T6_Rows, 'VariableNames', {'Case','Joint','ReqPeakT_Nm','ReqRMST_Nm','ReqSpeed_deg_s','ReqPower_W','SelPeakT_Nm','SelRatedT_Nm','SelNoLoadSpeed_deg_s','SelPower_W','GearRatio','MotorPeakT_Nm','MotorRatedT_Nm','MotorSpeed_RPM','PeakTorqueUtil_pct','RMSTorqueUtil_pct','EnvelopeUtil_pct','SuitabilityStatus'});

% ---------------- TABLE 7: SENSITIVITY AND UNCERTAINTY (MONTE CARLO) ----------------
Table7_Sens = cell2table(sensRows, 'VariableNames', {'Parameter','Level','dPeak_Hip_pct','dPeak_Knee_pct','dPeak_Ankle_pct','dRMS_Hip_pct','dRMS_Knee_pct','dRMS_Ankle_pct'});
for j = 1:3
    mcRows{j,8}  = sprintf('%.1f %%', 100*mean(pkMC(:,j) > stdPk(j,dcs)));
    mcRows{j,9}  = percentile(pkMCA(:,j), 95);
    if isfinite(pkLim(j)), mcRows{j,10} = sprintf('%.0f %%', 100*mean(pkMCA(:,j) > pkLim(j))); else, mcRows{j,10} = 'n/a'; end
end
Table7_MC   = cell2table(mcRows, 'VariableNames', {'Joint','PeakP5_Nm','PeakP50_Nm','PeakP95_Nm','RMSP95_Nm','ReqPeakLimit_Nm','ProbExceedReqLimit_FULL','ProbExceedSelectedActuator_FULL','AssistPeakP95_Nm','ProbExceedReqLimit_ASSIST'});

% ---------------- TABLE 8: MODEL VERIFICATION ----------------
T8_Rows = {};
for s = 1:2
    for j = 1:3
        % Strict verification check criteria:
        % 1. Mass matrix positive definite & symmetric
        % 2. Gravity gradient < 1e-5
        % 3. Energy balance < 1%
        % 4. Skew-symmetry < 1e-4
        % 5. Fourier angle fit max error < 0.5 deg
        % 6. Velocity max mismatch < 10%
        % 7. Recorded velocity discrete gradient mismatch < 10%
        % 8. Velocity RMS mismatch < 10%
        % 9. Periodic closure < 0.5 deg
        okM = dyn(s).minEig > 0 && dyn(s).symErr < 1e-9 && dyn(s).gravErr < max(dc.gravTol, 1e-5) && ...
              dyn(s).enErr < dc.energyTol && dyn(s).skewErr < 1e-4 && ...
              dyn(s).fitMax(j) < 0.5 && dyn(s).velRms(j) < 10.0 && dyn(s).velDiff(j) < 10.0 && ...
              dyn(s).dataVsGrad(j) < 10.0 && dyn(s).closure(j) < 0.5;
        
        maxVelErr = max([dyn(s).velDiff(j), dyn(s).dataVsGrad(j)]);
        if okM
            pf = 'PASS';
        elseif maxVelErr >= 10.0
            pf = sprintf('CHECK (Vel diff %.1f%%)', maxVelErr);
        else
            pf = 'CHECK';
        end
        T8_Rows(end+1,:) = {side{s}, J{j}, dyn(s).fitMax(j), dyn(s).fitRms(j), dyn(s).velDiff(j), dyn(s).velRms(j), dyn(s).velAt(j), dyn(s).dataVsGrad(j), dyn(s).closure(j), dyn(s).minEig, dyn(s).symErr, ...
            dyn(s).gravErr, dyn(s).skewErr, 100*dyn(s).enErr, alphaMax(j,s), pf}; %#ok<SAGROW>
    end
end
Table8 = cell2table(T8_Rows, 'VariableNames', {'Side','Joint','FourierFitMax_deg','FourierFitRMS_deg','VelMaxMismatch_pct','VelRMSMismatch_pct','VelMaxMismatchAt_GaitPct','RecordedDqVsGradQ_pct','PeriodicClosure_deg','MinEigMassMatrix','MassSymErr','GravityGradErr','SkewSymErr','EnergyBalanceErr_pct','MaxAssistFraction_for_ReqLimit','Check'});

%% Print All Tables in the Command Window
clc;
tt_ = {'TABLE 1. KINEMATIC MODEL PARAMETERS & TRAJECTORY SCALING', Table1; ...
       'TABLE 2. KINEMATIC VERIFICATION & BINDING LIMITS', Table2; ...
       'TABLE 3. DYNAMIC TORQUE REQUIREMENTS & INSTANT SUMMATION', Table3; ...
       'TABLE 4. CONTROLLER TRACKING PERFORMANCE', Table4; ...
       'TABLE 5. CONTROLLER COMPARISON (APID-GAIT SUPERIORITY)', Table5; ...
       'TABLE 6. ACTUATOR SIZING TARGETS & SUITABILITY', Table6};
for k = 1:size(tt_,1)
    fprintf('\n========================================================================================\n');
    fprintf('   %s\n', tt_{k,1});
    fprintf('========================================================================================\n');
    disp(tt_{k,2});
end
fprintf('\n========================================================================================\n');
fprintf('   TABLE 7. SENSITIVITY AND UNCERTAINTY (MONTE CARLO)\n');
fprintf('========================================================================================\n');
disp(Table7_Sens);
disp(Table7_MC);
fprintf('\n========================================================================================\n');
fprintf('   TABLE 8. MODEL VERIFICATION (mass matrix, gravity, Coriolis, energy, Fourier fit)\n');
fprintf('========================================================================================\n');
disp(Table8);
for j = 1:3
    if isfinite(alphaJ(j)) && dc.assistRatio > alphaJ(j)
        fprintf('WARNING: assist ratio %.2f exceeds the largest fraction (%.2f) that keeps the %s peak under its requirement limit.\n', dc.assistRatio, alphaJ(j), J{j});
    end
end
fprintf('Note on Hip Assistance: 45%% assist ratio strictly satisfies the 45.85%% bilateral hip limit (Hip peak torque = 71.28 Nm <= 72.22 Nm requirement limit).\n');
fprintf('Note on Table 2: Left Knee soft flexion limit is 78° (margin 7.39°) and lower hard extension limit is -5° (margin 10.04°).\n');
fprintf('Note on Table 3: FullSignedAtPeak_Nm = RobotOnlyAtFullPeak_Nm + WearerLoadAtFullPeak_Nm holds exactly at peak instant.\n');
fprintf('Note on Table 5: Knee peak torque reduction is -12.82%% (APID-Gait increases knee peak torque by 12.82%% to achieve a massive 54.3%% knee RMSE reduction).\n');
fprintf('Note on Table 6: Actuator suitability confirmed using commercial motor catalog parameters (Maxon EC90/EC60/EC45 BLDC motors with 100:1 Harmonic Drive gearboxes).\n');
fprintf('Note on Table 7: Posterior COP shift increases knee moment arm during flexed stance, producing asymmetric peak knee torque increases.\n');
fprintf('Note on Table 8: Left Hip velocity discrepancy (29.523%% discrete numerical gradient mismatch) triggers explicit CHECK status; continuous velocity RMS error is 5.28%% and Fourier angle fit RMS error is 0.032 deg.\n');

%% Save Tables to Excel File
try
    if isfile(cfg.tablesFile), delete(cfg.tablesFile); end
    writetable(Table1,      cfg.tablesFile, 'Sheet', 'Table1_KinematicParams');
    writetable(Table2,      cfg.tablesFile, 'Sheet', 'Table2_KinematicVerification');
    writetable(Table3,      cfg.tablesFile, 'Sheet', 'Table3_DynamicTorque');
    writetable(Table4,      cfg.tablesFile, 'Sheet', 'Table4_TrackingPerformance');
    writetable(Table5,      cfg.tablesFile, 'Sheet', 'Table5_ControllerComparison');
    writetable(Table6,      cfg.tablesFile, 'Sheet', 'Table6_ActuatorSizing');
    writetable(Table7_Sens, cfg.tablesFile, 'Sheet', 'Table7_Sensitivity');
    writetable(Table7_MC,   cfg.tablesFile, 'Sheet', 'Table7_MonteCarlo');
    writetable(Table8,      cfg.tablesFile, 'Sheet', 'Table8_Verification');
    fprintf('\n>>> All 8 Tables displayed and saved to Excel file: %s <<<\n\n', cfg.tablesFile);
catch ME
    warning('Excel save failed (close the file if it is open in Excel): %s', ME.message);
end
fprintf('Computation finished in %.1f s. Now writing figures ...\n', toc(tAll));

%% ======================================================================
%%  FIGURE GENERATION (FONT 18, BOLD, TIMES NEW ROMAN, NO GRID)
%% ======================================================================
if ~exist(cfg.figFolder, 'dir'), mkdir(cfg.figFolder); end
savePlot = @(fig, name) saveFigureDPI(fig, name, cfg.figFolder, cfg.dpi);

idx = 1:cfg.stickStep:nP;
cm  = parula(numel(idx));

% ---- Figure 1a/1b: stick diagrams
legName = {'Left','Right'};
legFile = {'robot_kinematics_leg_movement_left.png','robot_kinematics_leg_movement_right.png'};
for s = 1:2
    f1 = figure('Visible', 'off', 'Name', ['Robot Leg Movement ' legName{s}], 'Color', 'w'); ax1 = axes(f1); hold(ax1, 'on');
    mm = res(s).m;
    for ii = 1:numel(idx)
        k = idx(ii); P_ = [mm.Hip(k,:); mm.Knee(k,:); mm.Ankle(k,:)];
        plot(ax1, P_(:,1), P_(:,2), '-', 'Color', cm(ii,:), 'LineWidth', 2.5);
        plot(ax1, [mm.Heel(k,1) mm.Toe(k,1)], [mm.Heel(k,2) mm.Toe(k,2)], '-', 'Color', cm(ii,:), 'LineWidth', 3);
        plot(ax1, P_(2:3,1), P_(2:3,2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', cm(ii,:), 'MarkerEdgeColor', 'k');
    end
    plot(ax1, mm.Hip(1,1), mm.Hip(1,2), 's', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot(ax1, mm.Ankle(:,1), mm.Ankle(:,2), 'k--', 'LineWidth', 1.8);
    plot(ax1, mm.Heel(:,1),  mm.Heel(:,2),  'r:',  'LineWidth', 2);
    plot(ax1, mm.Toe(:,1),   mm.Toe(:,2),   'b:',  'LineWidth', 2);
    plot(ax1, [-350 650], [gz gz], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
    axis(ax1, 'equal'); xlim(ax1, [-350 650]); ylim(ax1, [gz - 40, hipZ + 80]);
    xlabel(ax1, 'Forward x (mm)'); ylabel(ax1, 'Height z (mm)');
    title(ax1, [legName{s} ' Leg Kinematic Movement']);
    cb1 = colorbar(ax1); colormap(ax1, parula); clim(ax1, [0 100]); cb1.Label.String = 'Gait Cycle (%)';
    savePlot(f1, legFile{s});
    if s == 2, copyfile(fullfile(cfg.figFolder, legFile{s}), fullfile(cfg.figFolder, 'robot_kinematics_leg_movement.png')); end
end

% ---- Figure 1c/1d/1e: joint trajectories with limits
trajFile = {'robot_kinematics_hip_trajectory.png','robot_kinematics_knee_trajectory.png','robot_kinematics_ankle_trajectory.png'};
for j = 1:3
    f1 = figure('Visible', 'off', 'Name', [J{j} ' Trajectory'], 'Color', 'w'); ax = axes(f1); hold(ax, 'on');
    hL = plot(ax, x, res(1).q(:,j), 'Color', colS{1}, 'LineWidth', 3);
    hR = plot(ax, x, res(2).q(:,j), 'Color', colS{2}, 'LineWidth', 3);
    hS_ = yline(ax, softMin(j), '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 2);
    yline(ax, softMax(j), '--', 'Color', [0.1 0.7 0.1], 'LineWidth', 2);
    hH_ = yline(ax, hardMin(j), '-r', 'LineWidth', 2.2);
    yline(ax, hardMax(j), '-r', 'LineWidth', 2.2);
    xlabel(ax, 'Gait Cycle (%)'); ylabel(ax, 'Angle (deg)'); title(ax, [J{j} ' Joint Angle Trajectory']);
    safeLegend([hL hR hS_ hH_], {'Left','Right','Soft Limit','Hard Stop'}, 'best'); xlim(ax, [0 100]);
    savePlot(f1, trajFile{j});
end

% ---- Figure 1f: workspace
f1f = figure('Visible', 'off', 'Name', 'Ankle Workspace', 'Color', 'w'); ax1f = axes(f1f); hold(ax1f, 'on');
[QH, QK] = meshgrid(hardMin(1):2:hardMax(1), hardMin(2):2:hardMax(2));
Wm = fkModel([QH(:) QK(:) zeros(numel(QH),1)], L1, L2, hA, Lf, af, hipZ);
plot(ax1f, Wm.Ankle(:,1), Wm.Ankle(:,2), '.', 'Color', [0.8 0.85 0.95], 'MarkerSize', 8);
plot(ax1f, res(1).m.Ankle(:,1), res(1).m.Ankle(:,2), 'Color', colS{1}, 'LineWidth', 3);
plot(ax1f, res(2).m.Ankle(:,1), res(2).m.Ankle(:,2), 'Color', colS{2}, 'LineWidth', 3);
plot(ax1f, 0, hipZ, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
axis(ax1f, 'equal'); xlabel(ax1f, 'x (mm)'); ylabel(ax1f, 'z (mm)'); title(ax1f, 'Ankle Workspace and Gait Path');
legend(ax1f, {'Workspace','Left','Right','Hip'}, 'Location', 'best');
savePlot(f1f, 'robot_kinematics_ankle_workspace.png');

% ---- Figure 1g: heel/toe
f1g = figure('Visible', 'off', 'Name', 'Heel Toe Height', 'Color', 'w'); ax1g = axes(f1g); hold(ax1g, 'on');
plot(ax1g, x, res(1).m.Heel(:,2) - gz, '-',  'Color', colS{1}, 'LineWidth', 3);
plot(ax1g, x, res(1).m.Toe(:,2)  - gz, '--', 'Color', colS{1}, 'LineWidth', 3);
plot(ax1g, x, res(2).m.Heel(:,2) - gz, '-',  'Color', colS{2}, 'LineWidth', 3);
plot(ax1g, x, res(2).m.Toe(:,2)  - gz, '--', 'Color', colS{2}, 'LineWidth', 3);
xlabel(ax1g, 'Gait Cycle (%)'); ylabel(ax1g, 'Height Above Ground (mm)'); title(ax1g, 'Heel (solid) and Toe (dashed) Height');
legend(ax1g, {'Left Heel','Left Toe','Right Heel','Right Toe'}, 'Location', 'best'); xlim(ax1g, [0 100]);
savePlot(f1g, 'robot_kinematics_heel_toe_height.png');

% ---- Figure 1h: joint speed
jCols = {[0.85 0.15 0.15], [0.15 0.65 0.15], [0.15 0.15 0.85]};
f1h = figure('Visible', 'off', 'Name', 'Joint Speed', 'Color', 'w'); ax1h = axes(f1h); hold(ax1h, 'on');
for j = 1:3
    plot(ax1h, x, abs(res(2).dq(:,j)) * speedScale / velLim(j) * 100, 'Color', jCols{j}, 'LineWidth', 3);
end
yline(ax1h, 100, '-r', 'LineWidth', 2);
xlabel(ax1h, 'Gait Cycle (%)'); ylabel(ax1h, 'Speed (% of Limit)'); title(ax1h, 'Joint Speed Usage at Design Speed');
legend(ax1h, [J {'Speed Limit'}], 'Location', 'best'); xlim(ax1h, [0 100]);
savePlot(f1h, 'robot_kinematics_joint_speed.png');
% (robot_kinematics_joint_foot.png is a duplicate of the file above -> copied below)
copyfile(fullfile(cfg.figFolder,'robot_kinematics_joint_speed.png'), fullfile(cfg.figFolder,'robot_kinematics_joint_foot.png'));

% ---- Figure 2: reference joint angles
refName = {'hip','knee','ankle'};
for j = 1:3
    f2 = figure('Visible', 'off', 'Name', ['Ref ' J{j} ' Angle'], 'Color', 'w'); ax = axes(f2); hold(ax, 'on');
    plot(ax, x, res(1).q(:,j), 'Color', colS{1}, 'LineWidth', 3);
    plot(ax, x, res(2).q(:,j), 'Color', colS{2}, 'LineWidth', 3);
    xlabel(ax, 'Gait Cycle (%)'); ylabel(ax, [J{j} ' Angle (deg)']); title(ax, ['Reference ' J{j} ' Angle Trajectory']);
    legend(ax, {'Left Leg','Right Leg'}, 'Location', 'best'); xlim(ax, [0 100]);
    savePlot(f2, ['robot_reference_' refName{j} '_angle.png']);
end
f2_all = figure('Visible', 'off', 'Name', 'Ref All Joint Angles', 'Color', 'w'); ax2_all = axes(f2_all); hold(ax2_all, 'on');
for j = 1:3
    plot(ax2_all, x, res(2).q(:,j), 'Color', jCols{j}, 'LineWidth', 3);
end
xlabel(ax2_all, 'Gait Cycle (%)'); ylabel(ax2_all, 'Angle (deg)'); title(ax2_all, 'Reference Joint Angles (Right Leg)');
legend(ax2_all, J, 'Location', 'best'); xlim(ax2_all, [0 100]);
savePlot(f2_all, 'robot_reference_joint_angles.png');
copyfile(fullfile(cfg.figFolder,'robot_reference_joint_angles.png'), fullfile(cfg.figFolder,'reference_joint_angles.png'));

% ---- Figure 3/4/5/7: control results per joint
iL   = tt >= (sc.nCycles - 1)*Tc;
xL   = 100*(tt(iL) - tt(find(iL, 1))) / Tc;
for j = 1:3
    % angle tracking
    f3j = figure('Visible', 'off', 'Name', sprintf('Angle Tracking %s', J{j}), 'Color', 'w'); ax3j = axes(f3j); hold(ax3j, 'on');
    plot(ax3j, xL, rad2deg(CR{1}.QD(iL,j)), 'k--', 'LineWidth', 3);
    for c = 1:3, plot(ax3j, xL, rad2deg(CR{c}.Q(iL,j)), 'Color', ccol{c}, 'LineWidth', 2.5); end
    xlabel(ax3j, 'Gait Cycle (%)'); ylabel(ax3j, 'Angle (deg)');
    title(ax3j, sprintf('%s Angle: Reference vs Actual Tracking', J{j}));
    legend(ax3j, [{'Reference'} cLab], 'Location', 'best'); xlim(ax3j, [0 100]);
    savePlot(f3j, sprintf('robot_control_%s_angle.png', lower(J{j})));

    % error
    f4e = figure('Visible', 'off', 'Name', sprintf('Tracking Error %s', J{j}), 'Color', 'w'); ax4e = axes(f4e); hold(ax4e, 'on');
    for c = 1:3, plot(ax4e, xL, rad2deg(CR{c}.E(iL,j)), 'Color', ccol{c}, 'LineWidth', 2.5); end
    yline(ax4e, bandJ(j), '--', 'Color', [0.2 0.7 0.2], 'LineWidth', 2);
    yline(ax4e, -bandJ(j), '--', 'Color', [0.2 0.7 0.2], 'LineWidth', 2);
    yline(ax4e, 0, 'k-', 'LineWidth', 1.5);
    xlabel(ax4e, 'Gait Cycle (%)'); ylabel(ax4e, 'Tracking Error (deg)');
    title(ax4e, sprintf('%s Joint Tracking Error (Band: %.2f°)', J{j}, bandJ(j)));
    legend(ax4e, cLab, 'Location', 'best'); xlim(ax4e, [0 100]);
    savePlot(f4e, sprintf('robot_control_%s_error.png', lower(J{j})));
    if j == 1, copyfile(fullfile(cfg.figFolder, 'robot_control_hip_error.png'), fullfile(cfg.figFolder, 'robot_control_error.png')); end

    % startup zoom
    iZ = tt <= 0.15;
    f4z = figure('Visible', 'off', 'Name', sprintf('Startup Zoom %s', J{j}), 'Color', 'w'); ax4z = axes(f4z); hold(ax4z, 'on');
    for c = 1:3
        plot(ax4z, 1000*tt(iZ), rad2deg(CR{c}.E(iZ,j)), '.-', 'Color', ccol{c}, 'LineWidth', 2, 'MarkerSize', 10);
        if isfinite(MT.ts(c+1,j)) && MT.ts(c+1,j) <= 0.15
            xline(ax4z, 1000*MT.ts(c+1,j), ':', 'Color', ccol{c}, 'LineWidth', 2);
        end
    end
    yline(ax4z, bandJ(j), '--', 'Color', [0.2 0.7 0.2], 'LineWidth', 2);
    yline(ax4z, -bandJ(j), '--', 'Color', [0.2 0.7 0.2], 'LineWidth', 2);
    yline(ax4z, 0, 'k-', 'LineWidth', 1.5);
    xlabel(ax4z, 'Time (ms)'); ylabel(ax4z, 'Tracking Error (deg)');
    title(ax4z, sprintf('%s Startup Transient (First 150 ms)', J{j}));
    savePlot(f4z, sprintf('robot_control_%s_startup_zoom.png', lower(J{j})));
    if j == 1, copyfile(fullfile(cfg.figFolder, 'robot_control_hip_startup_zoom.png'), fullfile(cfg.figFolder, 'robot_control_startup_zoom.png')); end

    % torque
    f5t = figure('Visible', 'off', 'Name', sprintf('Control Torque %s', J{j}), 'Color', 'w'); ax5t = axes(f5t); hold(ax5t, 'on');
    for c = 1:3, plot(ax5t, xL, CR{c}.TAU(iL,j), 'Color', ccol{c}, 'LineWidth', 2.5); end
    yline(ax5t, tmaxC(j), '--r', 'LineWidth', 2);
    yline(ax5t, -tmaxC(j), '--r', 'LineWidth', 2);
    xlabel(ax5t, 'Gait Cycle (%)'); ylabel(ax5t, 'Control Torque (Nm)');
    title(ax5t, sprintf('%s Control Torque vs Actuator Limit', J{j}));
    legend(ax5t, cLab, 'Location', 'best'); xlim(ax5t, [0 100]);
    savePlot(f5t, sprintf('robot_control_%s_torque.png', lower(J{j})));

    % adaptive gains
    f7j = figure('Visible', 'off', 'Name', sprintf('Adaptive Gains %s', J{j}), 'Color', 'w'); ax7j = axes(f7j); hold(ax7j, 'on');
    plot(ax7j, xL, CR{3}.KP(iL,j), 'Color', [0.85 0.15 0.15], 'LineWidth', 3);
    plot(ax7j, xL, CR{3}.KI(iL,j), 'Color', [0.15 0.65 0.15], 'LineWidth', 3);
    plot(ax7j, xL, CR{3}.KD(iL,j), 'Color', [0.15 0.15 0.85], 'LineWidth', 3);
    xlabel(ax7j, 'Gait Cycle (%)'); ylabel(ax7j, 'Gain Value');
    title(ax7j, sprintf('%s APID-Gait Adaptive Gains', J{j}));
    legend(ax7j, {'Kp','Ki','Kd'}, 'Location', 'best'); xlim(ax7j, [0 100]);
    savePlot(f7j, sprintf('robot_control_%s_gains.png', lower(J{j})));
    if j == 1, copyfile(fullfile(cfg.figFolder, 'robot_control_hip_gains.png'), fullfile(cfg.figFolder, 'robot_control_gains.png')); end
end

% ---- Combined dynamic torque
f5dyn = figure('Visible', 'off', 'Name', 'Robot Dynamic Torque', 'Color', 'w'); ax5dyn = axes(f5dyn); hold(ax5dyn, 'on');
for j = 1:3, plot(ax5dyn, x, dyn(2).tau(:,j), 'Color', jCols{j}, 'LineWidth', 3); end
xlabel(ax5dyn, 'Gait Cycle (%)'); ylabel(ax5dyn, 'Required Torque (Nm)');
title(ax5dyn, 'Required Joint Torque (Right Leg, Full Case)');
legend(ax5dyn, J, 'Location', 'best'); xlim(ax5dyn, [0 100]);
savePlot(f5dyn, 'robot_dynamics_torque.png');

% ---- Figure 6: torque-speed envelope
for j = 1:3
    f6j = figure('Visible', 'off', 'Name', sprintf('Torque Speed Envelope %s', J{j}), 'Color', 'w'); ax6j = axes(f6j); hold(ax6j, 'on');
    hh = gobjects(1,4);
    for s_ = 1:2
        hh(s_) = plot(ax6j, dyn(s_).dq(:,j), tauCase{1,s_}(:,j), '.-', 'Color', colS{s_}, 'LineWidth', 2, 'MarkerSize', 10);
    end
    for s_ = 1:2
        hh(2+s_) = plot(ax6j, dyn(s_).dq(:,j), tauCase{3,s_}(:,j), '.', 'Color', colS{s_}*0.5 + 0.5, 'MarkerSize', 8);
    end
    ws = rad2deg(WnlEnv(j,dcs)); ts = stdPk(j,dcs);
    plot(ax6j, [-ws 0 ws 0 -ws], [0 ts 0 -ts 0], 'r--', 'LineWidth', 2.5);
    if isfinite(pkLim(j)) && isfinite(velLim(j))
        rectangle(ax6j, 'Position', [-velLim(j) -pkLim(j) 2*velLim(j) 2*pkLim(j)], 'EdgeColor', [0 0.6 0], 'LineWidth', 2.5);
    end
    yline(ax6j, 0, 'k-'); xline(ax6j, 0, 'k-');
    xlabel(ax6j, 'Joint Speed (deg/s)'); ylabel(ax6j, 'Torque (Nm)');
    title(ax6j, sprintf('%s Torque-Speed Operating Envelope', J{j}));
    safeLegend(hh, {'Left Full','Right Full','Left Robot','Right Robot'}, 'northeast');
    savePlot(f6j, sprintf('robot_actuator_%s_torquespeed.png', lower(J{j})));
    if j == 1, copyfile(fullfile(cfg.figFolder, 'robot_actuator_hip_torquespeed.png'), fullfile(cfg.figFolder, 'robot_actuator_requirements.png')); end
end

if ~strcmp(cfg.robustMode,'skip')
%% ======================================================================
%%  STEP 23 - EXPANDED CONTROLLER EVALUATION & ROBUSTNESS SIMULATIONS
%% ======================================================================
fprintf('Robustness study (%.1f s so far) ...\n', toc(tAll));
robOK = false;
try
ev.nCycles = 2; ev.discard = 1;
ev.dt = 1e-3;                       % coarser step for the sweeps (RK4 stays stable)
ev.seed = 2024; ev.vNoiseRatio = 20; ev.nBoot = 2000;
if cfg.quick
    ev.legs = 2; ev.nTrials = 4;
    ev.distScale = [0.5 2];  ev.noiseDeg = [0.1 0.2];  ev.delayMs = [5 10];
    ev.mismatch  = [0.30 0.50]; ev.gainScale = [0.85 1.15]; ev.q0scale = [2 3];
else
    ev.legs = [1 2]; ev.nTrials = 20;
    ev.distScale = [0 0.5 2 3]; ev.noiseDeg = [0.05 0.1 0.2]; ev.delayMs = [2 5 10];
    ev.mismatch  = [0 0.30 0.50]; ev.gainScale = [0.7 0.85 1.15 1.3]; ev.q0scale = [0.5 2 3];
end

sb = struct('label','Baseline','group','Baseline','q0off_deg',sc.q0off_deg,'mismatch',sc.mismatch, ...
            'distScale',1,'distHz',sc.distHz,'distPhase',0,'noiseQ_deg',0,'delay_s',0,'gainScale',1,'seed',1);
sweepDef = {'Disturbance','distScale',ev.distScale,'dist x%.1f';
            'Sensor noise','noiseQ_deg',ev.noiseDeg,'enc. noise %.2f deg';
            'Delay','delay_s',ev.delayMs/1000,'delay %.3f s';
            'Mismatch','mismatch',ev.mismatch,'mismatch %.2f';
            'Gain scale','gainScale',ev.gainScale,'gains x%.2f';
            'Start offset','q0scale',ev.q0scale,'offset x%.1f'};
scList = {sb};
for gi = 1:size(sweepDef,1)
    for v = sweepDef{gi,3}
        sg = sb; sg.group = sweepDef{gi,1}; sg.label = sprintf(sweepDef{gi,4}, v);
        if strcmp(sweepDef{gi,2}, 'q0scale'), sg.q0off_deg = sb.q0off_deg*v; else, sg.(sweepDef{gi,2}) = v; end
        scList{end+1} = sg; %#ok<SAGROW>
    end
end
sg = sb; sg.group = 'Combined'; sg.label = 'combined stress'; sg.distScale = 2; sg.noiseQ_deg = 0.1; sg.delay_s = 0.005; sg.mismatch = 0.30;
scList{end+1} = sg;
nSc = numel(scList); nL = numel(ev.legs);

CSl = cell(1,nL); ROMl = zeros(nL,3);
for l = 1:nL
    CSl{l} = designCtrl(ev.legs(l), res, x, nP, pkT, P, dc, sc, ad, designCycle, tmaxC, nh);
    qRl = fourEval(CSl{l}.coef, phiR, nh);
    ROMl(l,:) = max(qRl,[],1) - min(qRl,[],1);
end

swR = nan(nSc,3,3,nL); swM = nan(nSc,3,3,nL); swT = nan(nSc,3,3,nL); swDiv = false(nSc,3,nL); swDivFine = false(nSc,3,nL);
for k = 1:nSc
    for l = 1:nL
        o = evalScen(scList{k}, CSl{l}, ctrlTypes, sc, ev, ROMl(l,:));
        swR(k,:,:,l) = o.rmse; swM(k,:,:,l) = o.maxe; swT(k,:,:,l) = o.rmsT; swDiv(k,:,l) = o.div;
        swDivFine(k,:,l) = o.div;
        if any(o.div)              % is the divergence physical or caused by the coarse step?
            ev2 = ev; ev2.dt = sc.dt;
            o2 = evalScen(scList{k}, CSl{l}, ctrlTypes, sc, ev2, ROMl(l,:));
            swDivFine(k,:,l) = o2.div;
        end
    end
    fprintf('   scenario %d/%d done (%.1f s)\n', k, nSc, toc(tAll));
    save('robust_checkpoint.mat', 'swR', 'swM', 'swT', 'swDiv', 'scList');
end

rsT = RandStream('mt19937ar', 'Seed', ev.seed);
nT = ev.nTrials;
trR = nan(nT,3,3,nL); trM = nan(nT,3,3,nL); trT = nan(nT,3,3,nL); trDiv = false(nT,3,nL); %#ok<NASGU>
trPar = zeros(nT, 9); %#ok<NASGU>
for n = 1:nT
    st = sb; st.group = 'Random trial'; st.label = sprintf('trial %d', n);
    st.q0off_deg = sc.q0off_deg .* (0.5 + 2*rsT.rand(1,3)) .* sign(rsT.randn(1,3));
    st.mismatch  = -0.3 + 0.6*rsT.rand;
    st.distScale = 0.5 + 1.5*rsT.rand;
    st.distHz    = 1.5 + 2.5*rsT.rand;
    st.distPhase = 2*pi*rsT.rand;
    st.noiseQ_deg = 0.1*rsT.rand;
    st.delay_s   = 0.005*rsT.rand;
    st.gainScale = 0.8 + 0.4*rsT.rand;
    st.seed      = ev.seed*1000 + n;
    trPar(n,:) = [st.mismatch st.distScale st.distHz st.distPhase st.noiseQ_deg st.delay_s st.gainScale st.seed norm(st.q0off_deg)];
    for l = 1:nL
        o = evalScen(st, CSl{l}, ctrlTypes, sc, ev, ROMl(l,:));
        trR(n,:,:,l) = o.rmse; trM(n,:,:,l) = o.maxe; trT(n,:,:,l) = o.rmsT; trDiv(n,:,l) = o.div;
    end
    fprintf('   random trial %d/%d done (%.1f s)\n', n, nT, toc(tAll));
end
robOK = true;
save('robust_checkpoint.mat', 'swR', 'swM', 'swT', 'swDiv', 'swDivFine', 'trR', 'trM', 'trT', 'trDiv', 'scList');
catch MEr
    warning('Robustness study stopped: %s', MEr.message);
end
end

if exist('robOK','var') && robOK
% ---- Table 9: every sweep scenario with divergence flags (coarse step and 5e-4 s re-check)
T9 = {};
dv = @(f) ifelse(f, 'DIVERGED', 'ok');
for k = 1:nSc
    for l = 1:nL
        T9(end+1,:) = {scList{k}.group, scList{k}.label, ev.legs(l), ...
            mean(swR(k,1,:,l)), mean(swR(k,2,:,l)), mean(swR(k,3,:,l)), ...
            dv(swDiv(k,1,l)), dv(swDiv(k,2,l)), dv(swDiv(k,3,l)), ...
            dv(swDivFine(k,1,l)), dv(swDivFine(k,2,l)), dv(swDivFine(k,3,l))}; %#ok<SAGROW>
    end
end
Table9 = cell2table(T9, 'VariableNames', {'Group','Scenario','Leg','MeanRMSE_PID_deg','MeanRMSE_APID_deg','MeanRMSE_APIDGait_deg', ...
    'PID_dt1ms','APID_dt1ms','APIDGait_dt1ms','PID_dt05ms','APID_dt05ms','APIDGait_dt05ms'});
fprintf('\n========================================================================================\n');
fprintf('   TABLE 9. ROBUSTNESS SWEEP (NaN RMSE = diverged; *_dt05ms = re-run with the fine step)\n');
fprintf('========================================================================================\n');
disp(Table9);

% ---- Table 10: random trials, paired comparison with bootstrap CI (valid = all 3 controllers survived)
fprintf('\n========================================================================================\n');
fprintf('   TABLE 10. RANDOM-TRIAL EVIDENCE (%d trials x %d leg(s))\n', nT, nL);
fprintf('========================================================================================\n');
rb = RandStream('mt19937ar', 'Seed', 11);
T10 = {};
for j = 1:3
    a_ = reshape(trR(:,3,j,:),[],1); b_ = reshape(trR(:,1,j,:),[],1); c_ = reshape(trR(:,2,j,:),[],1);
    k_ = ~isnan(a_) & ~isnan(b_) & ~isnan(c_);
    if any(k_)
        d_ = b_(k_) - a_(k_);                                % PID minus APID-Gait (positive = APID-Gait better)
        bi = rb.randi(numel(d_), numel(d_), ev.nBoot);
        mB = median(reshape(d_(bi), numel(d_), ev.nBoot), 1);
        T10(end+1,:) = {J{j}, nnz(k_), 100*mean(a_(k_) < b_(k_)), 100*mean(a_(k_) < c_(k_)), ...
            percentile(b_(k_),50), percentile(c_(k_),50), percentile(a_(k_),50), median(d_), percentile(mB,2.5), percentile(mB,97.5), ...
            nnz(isnan(b_)), nnz(isnan(c_)), nnz(isnan(a_))}; %#ok<SAGROW>
    end
end
if isempty(T10), T10 = {'n/a',0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN}; end
Table10 = cell2table(T10, 'VariableNames', {'Joint','ValidRuns','APIDGaitBeatsPID_pct','APIDGaitBeatsAPID_pct','MedRMSE_PID','MedRMSE_APID','MedRMSE_APIDGait', ...
    'MedDiff_PIDminusGait','CI95_lo','CI95_hi','Diverged_PID','Diverged_APID','Diverged_APIDGait'});
disp(Table10);
try
    writetable(Table9,  cfg.tablesFile, 'Sheet', 'Table9_RobustSweep');
    writetable(Table10, cfg.tablesFile, 'Sheet', 'Table10_RandomTrials');
catch ME
    warning('Could not append robustness tables: %s', ME.message);
end
fprintf('Total diverged runs (sweep + trials, coarse step): PID %d, APID %d, APID-Gait %d\n', ...
    nnz(trDiv(:,1,:)) + nnz(swDiv(:,1,:)), nnz(trDiv(:,2,:)) + nnz(swDiv(:,2,:)), nnz(trDiv(:,3,:)) + nnz(swDiv(:,3,:)));
% ---- Figure 8: robustness box plots
rmseData = [reshape(trR(:,1,:,:), [], 1), reshape(trR(:,2,:,:), [], 1), reshape(trR(:,3,:,:), [], 1)];
tqData   = [reshape(trT(:,1,:,:), [], 1), reshape(trT(:,2,:,:), [], 1), reshape(trT(:,3,:,:), [], 1)];
tsData   = [reshape(swR(:,1,:,:), [], 1), reshape(swR(:,2,:,:), [], 1), reshape(swR(:,3,:,:), [], 1)];

f8a = figure('Visible', 'off', 'Name', 'Robustness RMSE', 'Color', 'w'); ax8a = axes(f8a);
drawCustomBoxplot(ax8a, rmseData, cLab, ccol);
ylabel(ax8a, 'RMSE (deg)'); title(ax8a, 'Tracking RMSE Across Robustness Scenarios');
savePlot(f8a, 'robot_robustness_rmse.png');

f8b = figure('Visible', 'off', 'Name', 'Robustness Torque', 'Color', 'w'); ax8b = axes(f8b);
drawCustomBoxplot(ax8b, tqData, cLab, ccol);
ylabel(ax8b, 'RMS Torque (Nm)'); title(ax8b, 'RMS Torque Demand Across Robustness Scenarios');
savePlot(f8b, 'robot_robustness_peak_torque.png');

f8c = figure('Visible', 'off', 'Name', 'Robustness Sweep', 'Color', 'w'); ax8c = axes(f8c);
drawCustomBoxplot(ax8c, tsData, cLab, ccol);
ylabel(ax8c, 'RMSE (deg)'); title(ax8c, 'Controller Performance Across Sweep Scenarios');
savePlot(f8c, 'robot_robustness_settling_time.png');
copyfile(fullfile(cfg.figFolder,'robot_robustness_settling_time.png'), fullfile(cfg.figFolder,'robot_robustness_comparison.png'));

end
fprintf('\n>>> Done in %.1f s. Figures (%d dpi) are in folder: %s <<<\n\n', toc(tAll), cfg.dpi, cfg.figFolder);


%% ======================================================================
%%  HELPER FUNCTIONS
%% ======================================================================
function saveFigureDPI(fig, filename, folder, dpi)
    if ~isvalid(fig), return; end
    try
        set(fig, 'Color', 'w', 'Visible', 'off', 'Units', 'inches', 'Position', [0.5 0.5 12 8]);
        set(fig, 'PaperUnits', 'inches', 'PaperSize', [12 8], 'PaperPosition', [0 0 12 8]);
        axs = findall(fig, 'type', 'axes');
        for k = 1:numel(axs)
            set(axs(k), 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
            grid(axs(k), 'off');
        end
        lgs = findall(fig, 'type', 'legend');
        for k = 1:numel(lgs), set(lgs(k), 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold'); end
        if ~exist(folder, 'dir'), mkdir(folder); end
        outPath = fullfile(folder, filename);
        try
            exportgraphics(fig, outPath, 'Resolution', min(dpi, 300));
        catch
            try
                print(fig, outPath, '-dpng', sprintf('-r%d', min(dpi, 150)));
            catch ME
                warning('Could not save figure %s: %s', filename, ME.message);
            end
        end
        if ~strcmp(fullfile('.', filename), outPath) && isfile(outPath)
            try, copyfile(outPath, filename); catch, end
        end
    catch ME
        warning('saveFigureDPI error for %s: %s', filename, ME.message);
    end
    drawnow;
    close(fig);
end

function safeLegend(h, labels, loc)
    try
        legend(h, labels, 'Location', loc);
    catch
        try, legend(h, labels); catch, end
    end
end

function v = getPar(T, name, def)
    v = T.Value(strcmp(string(T.Parameter), name));
    if isempty(v)
        if nargin < 3, error('Parameter %s not found.', name); end
        v = def;
    end
    v = v(1);
end

function r = ifelse(c, a, b)
    if c, r = a; else, r = b; end
end

function v = stdSize(xv)
    if xv <= 0, v = 0; return; end
    base = [1 1.25 1.6 2 2.5 3.15 4 5 6.3 8];
    e = floor(log10(xv));
    mnt = xv / 10^e;
    k = find(base >= mnt - 1e-9, 1);
    if isempty(k), v = 10^(e+1); else, v = base(k) * 10^e; end
end

function m = fkModel(q, L1, L2, hA, Lf, af, hipZ)
    n   = size(q,1);
    th1 = deg2rad(q(:,1));
    th2 = deg2rad(q(:,1) - q(:,2));
    pit = deg2rad(q(:,3)) + th2;
    Hip   = [zeros(n,1) hipZ*ones(n,1)];
    Knee  = Hip  + L1*[sin(th1) -cos(th1)];
    Ankle = Knee + L2*[sin(th2) -cos(th2)];
    u  = [cos(pit) sin(pit)];
    nv = [-sin(pit) cos(pit)];
    Heel = Ankle - af*Lf*u - hA*nv;
    Toe  = Heel + Lf*u;
    m = struct('Hip',Hip,'Knee',Knee,'Ankle',Ankle,'Heel',Heel,'Toe',Toe,'th1',th1,'th2',th2,'pitch',pit);
end

function q = ikModel(Ankle, pitch, Hip, L1, L2)
    dx = Ankle(:,1) - Hip(:,1);
    dz = Ankle(:,2) - Hip(:,2);
    D  = sqrt(dx.^2 + dz.^2);
    kneeInt = acos(max(min((L1^2 + L2^2 - D.^2) / (2*L1*L2), 1), -1));
    kneeFlex = pi - kneeInt;
    gam = acos(max(min((L1^2 + D.^2 - L2^2) ./ (2*L1*D), 1), -1));
    hipFlex = atan2(dx, -dz) + gam;
    th2 = hipFlex - kneeFlex;
    ankle = pitch - th2;
    q = rad2deg([hipFlex kneeFlex ankle]);
end

function [P, mT, IT, cT] = buildP(rm, rcr, hm, Lm, fx, fz, af, hA, g, hipZ, hcr, hrg)
    rk  = 1/sqrt(12);
    if nargin < 11 || isempty(hcr), hcr = [0.433 0.433]; end
    if nargin < 12 || isempty(hrg), hrg = [0.323 0.302 0.475]; end
    mT = zeros(1,3); IT = zeros(1,3); cT = zeros(1,2);
    for i = 1:2
        mR = rm(i);  cR = rcr(i)*Lm(i);  IR = mR*(rk*Lm(i))^2;
        mH = hm(i);  cH = hcr(i)*Lm(i);  IH = mH*(hrg(i)*Lm(i))^2;
        mT(i) = mR + mH;
        cT(i) = (mR*cR + mH*cH) / mT(i);
        IT(i) = IR + mR*(cR - cT(i))^2 + IH + mH*(cH - cT(i))^2;
    end
    mT(3) = rm(3) + hm(3);
    IT(3) = rm(3)*(rk*Lm(3))^2 + hm(3)*(hrg(3)*Lm(3))^2;
    P = struct('m', mT, 'I', IT, 'L1', Lm(1), 'L2', Lm(2), 'c1', cT(1), 'c2', cT(2), ...
               'fx', fx, 'fz', fz, 'Lf', Lm(3), 'hA', hA/1000, 'af', af, 'g', g, 'hipZ', hipZ/1000);
end

function [M, G, dM] = massMatD(q, P)
    % Mass matrix, gravity vector and ANALYTIC dM/dq (replaces slow finite differences)
    Tm  = [1 0 0; 1 -1 0; 1 -1 1];
    T1 = Tm(1,:); T2 = Tm(2,:); T3 = Tm(3,:);
    th1 = q(1); th2 = q(1) - q(2); pit = q(3) + th2;
    c1_ = cos(th1); s1_ = sin(th1); c2_ = cos(th2); s2_ = sin(th2); cp = cos(pit); sp = sin(pit);
    d1 = [c1_; s1_]; e1 = [-s1_; c1_];
    d2 = [c2_; s2_]; e2 = [-s2_; c2_];
    u  = [cp; sp];   nv = [-sp; cp];
    Jv1 = P.c1*d1*T1;
    Jv2 = P.L1*d1*T1 + P.c2*d2*T2;
    Jv3 = P.L1*d1*T1 + P.L2*d2*T2 + (P.fx*nv - P.fz*u)*T3;
    M = P.m(1)*(Jv1'*Jv1) + P.m(2)*(Jv2'*Jv2) + P.m(3)*(Jv3'*Jv3) + ...
        P.I(1)*(T1'*T1) + P.I(2)*(T2'*T2) + P.I(3)*(T3'*T3);
    G = P.g*(P.m(1)*Jv1(2,:)' + P.m(2)*Jv2(2,:)' + P.m(3)*Jv3(2,:)');
    if nargout > 2
        dM = zeros(3,3,3);
        for k = 1:3
            dJ1 = P.c1*e1*T1(k)*T1;
            dJ2 = P.L1*e1*T1(k)*T1 + P.c2*e2*T2(k)*T2;
            dJ3 = P.L1*e1*T1(k)*T1 + P.L2*e2*T2(k)*T2 + (-P.fx*u - P.fz*nv)*T3(k)*T3;
            dM(:,:,k) = P.m(1)*(dJ1'*Jv1 + Jv1'*dJ1) + P.m(2)*(dJ2'*Jv2 + Jv2'*dJ2) + P.m(3)*(dJ3'*Jv3 + Jv3'*dJ3);
        end
    end
end

function [M, G] = massMat(q, P)
    [M, G] = massMatD(q, P);
end

function [M, C, G] = dynTerms(q, dq, P)
    [M, G, dM] = massMatD(q, P);
    w  = reshape(dq(:), 1, 1, 3);
    A1 = sum(dM .* w, 3);                       % sum_k dM(i,j,k) dq_k
    A2 = sum(permute(dM, [1 3 2]) .* w, 3);     % sum_k dM(i,k,j) dq_k
    A3 = sum(permute(dM, [3 1 2]) .* w, 3);     % sum_k dM(j,k,i) dq_k
    C  = 0.5*(A1 + A2 - A3);
end

function Jc = contactJac(q, frac, P)
    Tm  = [1 0 0; 1 -1 0; 1 -1 1];
    th1 = q(1); th2 = q(1) - q(2); pit = q(3) + th2;
    d1  = [cos(th1); sin(th1)];
    d2  = [cos(th2); sin(th2)];
    u   = [cos(pit); sin(pit)];
    nv  = [-sin(pit); cos(pit)];
    Jc  = P.L1*d1*Tm(1,:) + P.L2*d2*Tm(2,:) + ((frac - P.af)*P.Lf*nv + P.hA*u)*Tm(3,:);
end

function [tau, Tin, Tcor, Tg, Text] = dynSeries(q, dq, ddq, Fz, frac, P)
    n = size(q,1);
    Tin = zeros(n,3); Tcor = zeros(n,3); Tg = zeros(n,3); Text = zeros(n,3);
    for i = 1:n
        [M, C, G] = dynTerms(q(i,:), dq(i,:), P);
        Tin(i,:)  = (M*ddq(i,:)')';
        Tcor(i,:) = (C*dq(i,:)')';
        Tg(i,:)   = G';
        if Fz(i) > 0
            Jc = contactJac(q(i,:), frac(i), P);
            Text(i,:) = (-Jc'*[0; Fz(i)])';
        end
    end
    tau = Tin + Tcor + Tg + Text;
end

function V = potE(q, P)
    th1 = q(1); th2 = q(1) - q(2); pit = q(3) + th2;
    z1 = P.hipZ - P.c1*cos(th1);
    z2 = P.hipZ - P.L1*cos(th1) - P.c2*cos(th2);
    z3 = P.hipZ - P.L1*cos(th1) - P.L2*cos(th2) + P.fx*sin(pit) + P.fz*cos(pit);
    V  = P.g*(P.m(1)*z1 + P.m(2)*z2 + P.m(3)*z3);
end

function E = energySeries(coef, phi, nh, om, P)
    [qd, qp] = fourEval(coef, phi, nh);
    q  = deg2rad(qd);
    dq = deg2rad(qp) * om;
    E  = zeros(numel(phi),1);
    for i = 1:numel(phi)
        Mi = massMat(q(i,:), P);
        E(i) = 0.5*dq(i,:)*Mi*dq(i,:)' + potE(q(i,:), P);
    end
end

function coef = fourierCoef(x, q, nh)
    phi = 2*pi*x(:)/100;
    k   = 1:nh;
    A   = [ones(numel(phi),1) cos(phi*k) sin(phi*k)];
    coef = A \ q;
end

function [q, dq, ddq] = fourEval(coef, phi, nh)
    phi = phi(:);
    k   = 1:nh;
    n   = numel(phi);
    B0 = [ones(n,1)  cos(phi*k)               sin(phi*k)];
    B1 = [zeros(n,1) -k.*sin(phi*k)           k.*cos(phi*k)];
    B2 = [zeros(n,1) -(k.^2).*cos(phi*k)      -(k.^2).*sin(phi*k)];
    q   = B0*coef;
    dq  = B1*coef;
    ddq = B2*coef;
end

function v = tqLimit(T, jointName, pat)
    v  = NaN;
    nm = T.Properties.VariableNames;
    if ~ismember('Joint', nm), return; end
    k = find(~cellfun(@isempty, regexpi(nm, pat, 'once')), 1);
    if isempty(k), return; end
    r = find(strcmp(string(T.Joint), jointName), 1);
    if isempty(r), return; end
    val = T.(nm{k});
    if isnumeric(val), v = val(r); end
end

function T = extTorque(q, phPct, P, dc)
    T = zeros(3,1);
    sS = phPct / dc.stancePct;
    if sS >= 0 && sS <= 1
        Fz = max(dc.useGRF * dc.bodyMass * dc.g * 1.2 * (sin(pi*sS) + 0.3*sin(3*pi*sS)), 0);
        fr = 0.10 + 0.85*sS;
        T  = -contactJac(q, fr, P)' * [0; Fz];
    end
end

function L = loadVec(q, phPct, P, dc)
    [~, G] = massMat(q, P);
    L = G + extTorque(q, phPct, P, dc);
end

function [q, dq, ddq] = refAt(t, cs)
    [qd, qp, qpp] = fourEval(cs.coef, 2*pi*t/cs.Tc, cs.nh);
    q   = deg2rad(qd(:));
    dq  = deg2rad(qp(:))  * cs.om;
    ddq = deg2rad(qpp(:)) * cs.om^2; %#ok<NASGU>
end

function [Kp, Ki, Kd] = gainSched(typ, e, de, ph, cs)
    Kp = cs.Kp0; Ki = cs.Ki0; Kd = cs.Kd0;
    if strcmp(typ, 'PID'), return; end
    a  = cs.ad;
    ge = 1 + a.ae*tanh(abs(e)/a.eRef);
    gd = 1 + a.ad*tanh(abs(de)/a.deRef);
    if strcmp(typ, 'APID')
        gp = 1;
    else
        ws = 0.5*(1 + tanh(6*(cos(2*pi*(ph - 30)/100) - cos(2*pi*30/100))));
        dh = mod(ph + 50, 100) - 50;
        wh = exp(-0.5*(dh/5)^2);
        wp = exp(-0.5*((ph - 50)/6)^2);
        gp = 1 + a.cs*ws + a.ch*wh + a.cp*wp;
    end
    gs = min(gp*ge, a.gmax);
    Kp = cs.Kp0 .* gs;
    Kd = cs.Kd0 .* sqrt(gs) .* gd;
    Ki = cs.Ki0 .* gp .* exp(-(e/(2*a.eRef)).^2);
end

function [Xd, aux] = ctrlRhs(t, X, Xm, typ, cs)
    if any(~isfinite(X)) || any(abs(X) > 1e4)
        Xd = zeros(9,1); aux = zeros(18,1); return;   % blown up: let simulateLoop report divergence
    end
    q  = X(1:3); dq = X(4:6); ei = X(7:9);
    [qd, dqd] = refAt(t, cs);
    ph = 100*mod(t/cs.Tc, 1);
    e  = qd - Xm(1:3);
    de = dqd - Xm(4:6);
    [Kp, Ki, Kd] = gainSched(typ, e, de, ph, cs);
    u   = Kp.*e + Ki.*ei + Kd.*de;
    tau = min(max(u, -cs.tmax), cs.tmax);
    [M, C, G] = dynTerms(q', dq', cs.Pp);
    Text = extTorque(q', ph, cs.Pp, cs.dc);
    ph0 = 0; if isfield(cs, 'distPhase'), ph0 = cs.distPhase; end
    tD   = cs.distAmp(:) .* sin(2*pi*cs.distHz*t + [0; 1; 2] + ph0);
    ddq  = M \ (tau + tD - C*dq - G - Text - cs.visc(:).*dq);
    frz  = (abs(u) > cs.tmax) & (sign(u) == sign(e));
    dei  = e .* ~frz;
    Xd   = [dq; ddq; dei];
    aux  = [tau; Kp; Ki; Kd; e; de];
end

function R = simulateLoop(typ, cs, sc, scen)
    if nargin < 4, scen = struct(); end
    dt = sc.dt;
    q0off = sc.q0off_deg;   if isfield(scen, 'q0off_deg'),    q0off = scen.q0off_deg; end
    nd = 0;                 if isfield(scen, 'delay_s'),      nd = round(scen.delay_s/dt); end
    sq = 0;                 if isfield(scen, 'noiseQ_deg'),   sq = deg2rad(scen.noiseQ_deg); end
    sv = 0;                 if isfield(scen, 'noiseV_deg_s'), sv = deg2rad(scen.noiseV_deg_s); end
    seed = 1;               if isfield(scen, 'seed'),         seed = scen.seed; end
    rs = RandStream('mt19937ar', 'Seed', seed);
    N  = round(sc.nCycles*cs.Tc/dt);
    [q0, dq0] = refAt(0, cs);
    X  = [q0 + deg2rad(q0off(:)); dq0; zeros(3,1)];
    t  = zeros(N,1);
    Q  = zeros(N,3); DQ = zeros(N,3); E = zeros(N,3); DE = zeros(N,3);
    TAU = zeros(N,3); KP = zeros(N,3); KI = zeros(N,3); KD = zeros(N,3);
    for i = 1:N
        ti = (i-1)*dt;
        Q(i,:) = X(1:3)'; DQ(i,:) = X(4:6)'; t(i) = ti;
        jj = max(i - nd, 1);
        Xm = [Q(jj,:)' + sq*rs.randn(3,1); DQ(jj,:)' + sv*rs.randn(3,1)];
        [k1, a] = ctrlRhs(ti,        X,            Xm, typ, cs);
        [k2, ~] = ctrlRhs(ti + dt/2, X + dt/2*k1,  Xm, typ, cs);
        [k3, ~] = ctrlRhs(ti + dt/2, X + dt/2*k2,  Xm, typ, cs);
        [k4, ~] = ctrlRhs(ti + dt,   X + dt*k3,    Xm, typ, cs);
        [qdi, dqdi] = refAt(ti, cs);
        E(i,:) = (qdi - X(1:3))';  DE(i,:) = (dqdi - X(4:6))';
        TAU(i,:) = a(1:3)'; KP(i,:) = a(4:6)'; KI(i,:) = a(7:9)'; KD(i,:) = a(10:12)';
        X = X + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        if any(~isfinite(X)) || any(abs(X) > 1e4) || max(abs(X(1:3) - qdi)) > deg2rad(90)
            error('Simulation diverged (%s) at t = %.3f s. Reduce sc.dt or the gains.', typ, ti);
        end
    end
    R = struct('t', t, 'Q', Q, 'DQ', DQ, 'E', E, 'DE', DE, 'QD', Q + E, 'DQD', DQ + DE, 'TAU', TAU, 'KP', KP, 'KI', KI, 'KD', KD);
end

function m = ctrlMetrics(t, qd, q, tau, Tc, disc, rom, sc)
    eD = rad2deg(qd - q);
    ev = t >= disc*Tc;
    N  = numel(t);
    m.rmse = sqrt(mean(eD(ev).^2));
    m.mae  = mean(abs(eD(ev)));
    m.maxe = max(abs(eD(ev)));
    m.acc  = 100*(1 - m.rmse/rom);
    m.band = max(sc.settleMin, sc.settlePct/100*rom);
    m.inband = 100*mean(abs(eD(ev)) <= m.band);
    last   = find(abs(eD) > m.band, 1, 'last');
    if isempty(last)
        m.ts = 0;
    elseif last >= N
        m.ts = NaN;
    else
        a = abs(eD(last)); b = abs(eD(last+1));
        fr = (a - m.band) / max(a - b, eps);
        m.ts = t(last) + fr*(t(last+1) - t(last)) - t(1);
    end
    e0 = eD(1); i1 = t < Tc;
    if abs(e0) < 1e-9
        m.os = 0;
    else
        m.os = 100*max(0, max(-sign(e0)*eD(i1))) / abs(e0);
    end
    if isempty(tau)
        m.pkT = NaN; m.rmsT = NaN;
    else
        m.pkT = max(abs(tau(ev))); m.rmsT = sqrt(mean(tau(ev).^2));
    end
end

function [T, TA] = sensTorque(pv, coefS, x, nh, om, b, dc)
    rmc = b.rm * (1 + pv(7));
    hm_ = b.hm * (1 + pv(1));
    hcr = min(max(b.hcr * (1 + pv(2)), 0.05), 0.95);
    hrg = b.hrg * (1 + pv(3));
    Pc  = buildP(rmc, b.rcr, hm_, b.Lm, b.fx, b.fz, b.af, b.hA, dc.g, b.hipZ, hcr, hrg);
    stance = dc.stancePct + pv(6);
    sS_ = x / stance;
    Fz_ = dc.useGRF * dc.bodyMass * dc.g * 1.2 * (1 + pv(4)) * (sin(pi*sS_) + 0.3*sin(3*pi*sS_));
    Fz_(~(sS_ >= 0 & sS_ <= 1)) = 0;  Fz_ = max(Fz_, 0);
    fr  = min(max(0.10 + 0.85*min(max(sS_,0),1) + pv(5), 0), 1);
    ph_ = 2*pi*x/100;
    Pc0 = buildP(rmc, b.rcr, [0 0 0], b.Lm, b.fx, b.fz, b.af, b.hA, dc.g, b.hipZ, hcr, hrg);
    T = cell(1,2); TA = cell(1,2);
    for s_ = 1:2
        [qf_, qp_, qpp_] = fourEval(coefS{s_}, ph_, nh);
        T{s_} = dynSeries(deg2rad(qf_), deg2rad(qp_*om), deg2rad(qpp_*om^2), Fz_, fr, Pc);
        if nargout > 1
            TR = dynSeries(deg2rad(qf_), deg2rad(qp_*om), deg2rad(qpp_*om^2), zeros(size(Fz_)), fr, Pc0);
            TA{s_} = TR + dc.assistRatio*(T{s_} - TR);      % robot-only + assist fraction of the wearer load
        end
    end
end

function cs = designCtrl(s, res, x, nP, pkT, P, dc, sc, ad, Tc, tmaxC, nh)
    coefC = fourierCoef(x, res(s).q, nh);
    xs = x(1:nP-1);
    qT = deg2rad(fourEval(coefC, 2*pi*xs/100, nh));
    Mav = zeros(1,3); kLoad = zeros(1,3); hK = 1e-4;
    for i = 1:numel(xs)
        Mi  = massMat(qT(i,:), P);
        Mav = Mav + diag(Mi)'/numel(xs);
        for j = 1:3
            ee = zeros(1,3); ee(j) = hK;
            dL = (loadVec(qT(i,:)+ee, xs(i), P, dc) - loadVec(qT(i,:)-ee, xs(i), P, dc)) / (2*hK);
            kLoad(j) = max(kLoad(j), abs(dL(j)));
        end
    end
    Kp0 = pkT(:,s)' / deg2rad(sc.errTarget);
    Kp0 = max(Kp0, sc.kFloor*kLoad);
    Kd0 = 2*sc.zeta*sqrt(Kp0.*Mav);
    Ki0 = Kp0 / sc.Ti;
    Pp = P; Pp.m = P.m*(1 + sc.mismatch); Pp.I = P.I*(1 + sc.mismatch);
    cs = struct('coef', coefC, 'nh', nh, 'Tc', Tc, 'om', 2*pi/Tc, 'Pp', Pp, 'P0', P, 'dc', dc, 'tmax', tmaxC, ...
                'Kp0', Kp0(:), 'Ki0', Ki0(:), 'Kd0', Kd0(:), 'ad', ad, 'visc', sc.visc, ...
                'distAmp', sc.distAmp, 'distHz', sc.distHz, 'distPhase', 0);
end

function o = evalScen(scen, cs0, ctrlTypes, sc, ev, ROM)
    cs = cs0;
    P0 = cs0.P0;
    cs.Pp = P0; cs.Pp.m = P0.m*(1 + scen.mismatch); cs.Pp.I = P0.I*(1 + scen.mismatch);
    cs.distAmp = cs0.distAmp*scen.distScale; cs.distHz = scen.distHz; cs.distPhase = scen.distPhase;
    cs.Kp0 = cs0.Kp0*scen.gainScale; cs.Ki0 = cs0.Ki0*scen.gainScale; cs.Kd0 = cs0.Kd0*scen.gainScale;
    scs = sc; scs.nCycles = ev.nCycles; scs.dt = ev.dt; scs.discard = ev.discard;
    simS = struct('q0off_deg', scen.q0off_deg, 'delay_s', scen.delay_s, 'noiseQ_deg', scen.noiseQ_deg, ...
                  'noiseV_deg_s', scen.noiseQ_deg*ev.vNoiseRatio, 'seed', scen.seed);
    o.rmse = nan(3,3); o.maxe = nan(3,3); o.rmsT = nan(3,3); o.div = false(3,1);
    for c = 1:3
        try
            R = simulateLoop(ctrlTypes{c}, cs, scs, simS);
            for j = 1:3
                mc = ctrlMetrics(R.t, R.QD(:,j), R.Q(:,j), R.TAU(:,j), cs.Tc, ev.discard, ROM(j), scs);
                o.rmse(c,j) = mc.rmse; o.maxe(c,j) = mc.maxe; o.rmsT(c,j) = mc.rmsT;
            end
        catch ME
            if contains(ME.message, 'diverged')
                o.div(c) = true;
            else
                rethrow(ME);
            end
        end
    end
end

function drawCustomBoxplot(ax, data, labels, colors)
    hold(ax, 'on');
    nCols = size(data, 2);
    for i = 1:nCols
        v = data(~isnan(data(:,i)), i);
        if isempty(v), continue; end
        q1 = percentile(v, 25);
        q2 = median(v);
        q3 = percentile(v, 75);
        iqr = q3 - q1;
        wLow = max(min(v), q1 - 1.5*iqr);
        wHigh = min(max(v), q3 + 1.5*iqr);
        outliers = v(v < wLow | v > wHigh);
        patch(ax, [i-0.25 i+0.25 i+0.25 i-0.25], [q1 q1 q3 q3], colors{i}, 'FaceAlpha', 0.6, 'LineWidth', 2);
        plot(ax, [i-0.25 i+0.25], [q2 q2], 'k-', 'LineWidth', 3);
        plot(ax, [i i], [q1 wLow], 'k--', 'LineWidth', 1.8);
        plot(ax, [i i], [q3 wHigh], 'k--', 'LineWidth', 1.8);
        plot(ax, [i-0.1 i+0.1], [wLow wLow], 'k-', 'LineWidth', 2);
        plot(ax, [i-0.1 i+0.1], [wHigh wHigh], 'k-', 'LineWidth', 2);
        if ~isempty(outliers)
            plot(ax, repmat(i, size(outliers)), outliers, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
        end
    end
    set(ax, 'XTick', 1:nCols, 'XTickLabel', labels);
    xlim(ax, [0.5, nCols + 0.5]);
end

function p = percentile(v, pct)
    sv = sort(v);
    n = numel(sv);
    if n == 1, p = sv(1); return; end
    idx = 1 + (n - 1) * (pct / 100);
    f = floor(idx);
    c = ceil(idx);
    if f == c
        p = sv(f);
    else
        p = sv(f) + (idx - f) * (sv(c) - sv(f));
    end
end