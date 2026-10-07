clc; clear; close all;

cfg.gaitFolder   = '.';
cfg.reqFile      = 'robot_requirements.xlsx';
cfg.designFile   = 'robot_design.xlsx';
cfg.outFile      = 'robot_kinematics.xlsx';
cfg.figPrefix    = 'robot_kinematics';
cfg.stickStep    = 5;
cfg.ikTol_deg    = 1e-6;
cfg.velTol       = 1.0;
cfg.swingWindow  = [70 90];

J    = {'Hip','Knee','Ankle'};
side = {'Left','Right'};
pre  = {'L','R'};
colS = {[0 0.447 0.741], [0.850 0.325 0.098]};

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

fprintf('\n========================================\n');
fprintf('   ROBOT KINEMATIC MODEL\n');
fprintf('========================================\n');
fprintf('Planar chain  : Hip -> Thigh (L1) -> Knee -> Shank (L2) -> Ankle -> Foot\n');
fprintf('Thigh L1      : %.2f mm (adjustable %.1f to %.1f)\n', L1, L1mn, L1mx);
fprintf('Shank L2      : %.2f mm (adjustable %.1f to %.1f)\n', L2, L2mn, L2mx);
fprintf('Ankle height  : %.2f mm | Foot plate : %.1f mm | Ankle at %.0f %% of foot from heel\n', hA, Lf, 100*af);
fprintf('Hip height    : %.1f mm (hip fixed in the model frame, x forward, z up)\n', hipZ);
fprintf('Cycle time    : recorded %.3f s | design %.3f s | speed scale x%.2f\n', recCycle, designCycle, speedScale);
fprintf('\nFORWARD KINEMATICS (flexion positive):\n');
fprintf('  th1   = qh              (thigh angle from vertical)\n');
fprintf('  th2   = qh - qk         (shank angle from vertical)\n');
fprintf('  pitch = qa + th2        (foot angle, toes up positive)\n');
fprintf('  Knee  = Hip   + L1*[ sin(th1), -cos(th1)]\n');
fprintf('  Ankle = Knee  + L2*[ sin(th2), -cos(th2)]\n');
fprintf('  Heel  = Ankle - af*Lf*u - hA*n,  Toe = Heel + Lf*u\n');
fprintf('Note: the knee term is qh - qk because the gait model defines knee flexion positive.\n');
fprintf('      The form qh + qk applies only if knee flexion is defined with the opposite sign.\n');

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

fprintf('\n========================================\n');
fprintf('   FORWARD KINEMATICS RESULTS (mm, deg)\n');
fprintf('========================================\n');
for s = 1:2
    m = res(s).m;
    fprintf('\n%s leg\n', side{s});
    fprintf('%-18s %10s %10s %10s\n', 'Quantity', 'Min', 'Max', 'Range');
    items = {'Knee x', m.Knee(:,1); 'Knee z', m.Knee(:,2); 'Ankle x', m.Ankle(:,1); 'Ankle z', m.Ankle(:,2); ...
             'Heel z (ground)', m.Heel(:,2) - gz; 'Toe z (ground)', m.Toe(:,2) - gz; 'Foot angle (deg)', rad2deg(m.pitch)};
    for k = 1:size(items,1)
        v = items{k,2};
        fprintf('%-18s %10.1f %10.1f %10.1f\n', items{k,1}, min(v), max(v), max(v) - min(v));
    end
    sw = x >= cfg.swingWindow(1) & x <= cfg.swingWindow(2);
    fprintf('Toe clearance in mid-swing (%d-%d %%) : %.1f mm\n', cfg.swingWindow(1), cfg.swingWindow(2), min(m.Toe(sw,2) - gz));
    fprintf('Ankle forward travel (hip frame)     : %.1f mm\n', max(m.Ankle(:,1)) - min(m.Ankle(:,1)));
end

fprintf('\n========================================\n');
fprintf('   JOINT LIMIT CHECK\n');
fprintf('========================================\n');
fprintf('%-6s %-6s %9s %9s %12s %12s %12s\n', 'Joint', 'Side', 'q min', 'q max', 'Soft margin', 'Hard margin', 'Outside soft');
softMargin = zeros(3,2); hardMargin = zeros(3,2); outSoft = zeros(3,2);
for s = 1:2
    for j = 1:3
        qj = res(s).q(:,j);
        softMargin(j,s) = min([qj - softMin(j); softMax(j) - qj]);
        hardMargin(j,s) = min([qj - hardMin(j); hardMax(j) - qj]);
        outSoft(j,s)    = 100 * mean(qj < softMin(j) | qj > softMax(j));
        fprintf('%-6s %-6s %9.2f %9.2f %11.2f° %11.2f° %11.1f %%\n', J{j}, side{s}, min(qj), max(qj), softMargin(j,s), hardMargin(j,s), outSoft(j,s));
    end
end

fprintf('\n========================================\n');
fprintf('   REACH, SINGULARITY AND JACOBIAN\n');
fprintf('========================================\n');
fprintf('%-6s %14s %14s %14s %14s %14s\n', 'Side', 'Max reach', 'Min knee flex', 'Min |det J|', 'Max cond(J)', 'Ankle v err');
vErr = zeros(1,2);
for s = 1:2
    vErr(s) = 100 * max(abs(res(s).spJ - res(s).spN)) / max(res(s).spN);
    fprintf('%-6s %13.2f%% %13.2f° %11.0f mm2 %14.1f %13.1f %%\n', side{s}, 100*max(res(s).reach), min(res(s).q(:,2)), ...
            min(abs(res(s).detJ)), max(res(s).condJ), vErr(s));
end
fprintf('det J = L1*L2*sin(qk): the leg is singular at full knee extension (qk = 0).\n');
fprintf('Ankle speed from J*dq vs finite difference of the ankle position: max difference shown above.\n');

fprintf('\n========================================\n');
fprintf('   FK / IK ROUND-TRIP AND LIMIT-CLIPPING TEST\n');
fprintf('========================================\n');
fprintf('%-6s %-6s %16s %26s\n', 'Side', 'Joint', 'Max FK->IK error', 'Ankle error if clipped to soft limit');
for s = 1:2
    for j = 1:3
        if j == 1
            fprintf('%-6s %-6s %13.2e ° %23.3f mm\n', side{s}, J{j}, max(res(s).ikErr(:,j)), max(res(s).clipErr));
        else
            fprintf('%-6s %-6s %13.2e °\n', side{s}, J{j}, max(res(s).ikErr(:,j)));
        end
    end
end

fprintf('\n========================================\n');
fprintf('   JOINT SPEED vs ACTUATOR REQUIREMENT\n');
fprintf('========================================\n');
fprintf('%-6s %-6s %14s %16s %14s %10s\n', 'Joint', 'Side', 'Peak (rec.)', 'Peak (design)', 'Limit', 'Usage');
velUse = zeros(3,2);
for s = 1:2
    for j = 1:3
        pk = max(abs(res(s).dq(:,j)));
        velUse(j,s) = pk * speedScale / velLim(j);
        fprintf('%-6s %-6s %10.1f d/s %12.1f d/s %10.1f d/s %9.0f %%\n', J{j}, side{s}, pk, pk*speedScale, velLim(j), 100*velUse(j,s));
    end
end

fprintf('\n========================================\n');
fprintf('   LINK-LENGTH ADJUSTMENT ANALYSIS (right leg)\n');
fprintf('========================================\n');
lens = [L1mn L2mn; L1 L2; L1mx L2mx];
lensName = {'Minimum', 'Nominal', 'Maximum'};
mNom = res(2).m;
Dnom = sqrt(sum((mNom.Ankle - mNom.Hip).^2, 2));
adj = struct();
fprintf('%-9s %8s %8s %22s %20s %18s\n', 'Setting', 'L1', 'L2', 'Ankle shift, same q', 'Reach of human path', 'Unreachable cycle');
for c = 1:3
    mc = fkModel(res(2).q, lens(c,1), lens(c,2), hA, Lf, af, hipZ);
    dev = sqrt(sum((mc.Ankle - mNom.Ankle).^2, 2));
    rr  = Dnom / (lens(c,1) + lens(c,2));
    adj(c).Ankle = mc.Ankle; adj(c).reach = rr;
    fprintf('%-9s %8.1f %8.1f %19.1f mm %18.1f %% %16.0f %%\n', lensName{c}, lens(c,1), lens(c,2), max(dev), 100*max(rr), 100*mean(rr > 1));
end
fprintf('Same joint angles on different link lengths move the ankle, so the links must be set to the wearer.\n');

allOK = all(outSoft(:) == 0) && all(cellfun(@(r) max(r(:)), {res.reach}) <= 1) && ...
        max(cellfun(@(e) max(e(:)), {res.ikErr})) < cfg.ikTol_deg && all(velUse(:) <= cfg.velTol);
fprintf('\n========================================\n');
fprintf('   VERIFICATION RESULT\n');
fprintf('========================================\n');
fprintf('Joint angles inside soft limits : %s\n', passStr(all(outSoft(:) == 0)));
fprintf('Gait path inside leg reach      : %s\n', passStr(all(cellfun(@(r) max(r(:)), {res.reach}) <= 1)));
fprintf('FK/IK round trip                : %s\n', passStr(max(cellfun(@(e) max(e(:)), {res.ikErr})) < cfg.ikTol_deg));
fprintf('Joint speed within requirement  : %s\n', passStr(all(velUse(:) <= cfg.velTol)));
if allOK
    fprintf('The robot can reproduce the extracted gait trajectory.\n');
else
    fprintf('At least one check failed. Review the tables above.\n');
end

% ---------------- KINEMATIC PLOTS ----------------
idx = 1:cfg.stickStep:nP;
cm  = parula(numel(idx));

f1 = figure('Name', 'Robot Leg Movement', 'Color', 'w', 'Position', [40 60 1500 650]);
for s = 1:2
    ax = subplot(1,2,s); hold(ax, 'on');
    m = res(s).m;
    for ii = 1:numel(idx)
        k = idx(ii);
        P = [m.Hip(k,:); m.Knee(k,:); m.Ankle(k,:)];
        plot(ax, P(:,1), P(:,2), '-', 'Color', cm(ii,:), 'LineWidth', 2.2);
        plot(ax, [m.Heel(k,1) m.Toe(k,1)], [m.Heel(k,2) m.Toe(k,2)], '-', 'Color', cm(ii,:), 'LineWidth', 3);
        plot(ax, P(2:3,1), P(2:3,2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', cm(ii,:), 'MarkerEdgeColor', 'k');
    end
    plot(ax, m.Hip(1,1), m.Hip(1,2), 's', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot(ax, m.Ankle(:,1), m.Ankle(:,2), 'k--', 'LineWidth', 1);
    plot(ax, m.Heel(:,1),  m.Heel(:,2),  'r:',  'LineWidth', 1.5);
    plot(ax, m.Toe(:,1),   m.Toe(:,2),   'b:',  'LineWidth', 1.5);
    plot(ax, [-350 650], [gz gz], '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
    axis(ax, 'equal'); grid(ax, 'on');
    xlim(ax, [-350 650]); ylim(ax, [gz - 40, hipZ + 80]);
    xlabel(ax, 'Forward x (mm)'); ylabel(ax, 'Height z (mm)');
    title(ax, sprintf('%s leg - robot configuration every %d %% of the gait cycle', side{s}, cfg.stickStep));
    colormap(ax, parula); caxis(ax, [0 100]);
    cb = colorbar(ax); cb.Label.String = 'Gait cycle (%)';
    if s == 1
        legend(ax, {'Leg links','Foot','Knee/ankle joints'}, 'Location', 'southwest');
    else
        text(ax, -330, hipZ + 50, 'Hip', 'FontWeight', 'bold');
        text(ax, 380, gz + 15, 'ankle path (dashed), heel (red), toe (blue)', 'FontSize', 8);
    end
end
sgtitle('Robot leg movement through one gait cycle');
saveas(f1, [cfg.figPrefix '_leg_movement.png']);

f2 = figure('Name', 'Joint and Foot Kinematics', 'Color', 'w', 'Position', [60 40 1500 800]);
for j = 1:3
    subplot(2,3,j); hold on;
    for s = 1:2
        plot(x, res(s).q(:,j), 'Color', colS{s}, 'LineWidth', 2);
    end
    yline(softMin(j), '--g', 'LineWidth', 1.2); yline(softMax(j), '--g', 'LineWidth', 1.2);
    yline(hardMin(j), '-r', 'LineWidth', 1.4);  yline(hardMax(j), '-r', 'LineWidth', 1.4);
    xlabel('Gait cycle (%)'); ylabel('deg'); title([J{j} ' angle']); grid on; xlim([0 100]);
    if j == 1, legend('Left','Right','Soft limit','','Hard stop','Location','best'); end
end
subplot(2,3,4); hold on;
[QH, QK] = meshgrid(hardMin(1):2:hardMax(1), hardMin(2):2:hardMax(2));
Wm = fkModel([QH(:) QK(:) zeros(numel(QH),1)], L1, L2, hA, Lf, af, hipZ);
plot(Wm.Ankle(:,1), Wm.Ankle(:,2), '.', 'Color', [0.82 0.86 0.95], 'MarkerSize', 6);
for s = 1:2
    plot(res(s).m.Ankle(:,1), res(s).m.Ankle(:,2), 'Color', colS{s}, 'LineWidth', 2);
end
plot(0, hipZ, 'ks', 'MarkerFaceColor', 'k');
axis equal; grid on; xlabel('x (mm)'); ylabel('z (mm)');
title('Ankle workspace and gait path');
legend('Workspace','Left','Right','Hip','Location','best');
subplot(2,3,5); hold on;
for s = 1:2
    plot(x, res(s).m.Heel(:,2) - gz, '-',  'Color', colS{s}, 'LineWidth', 2);
    plot(x, res(s).m.Toe(:,2) - gz,  '--', 'Color', colS{s}, 'LineWidth', 2);
end
xlabel('Gait cycle (%)'); ylabel('mm above ground'); title('Heel (solid) and toe (dashed) height'); grid on; xlim([0 100]);
subplot(2,3,6); hold on;
for j = 1:3
    plot(x, abs(res(2).dq(:,j)) * speedScale / velLim(j) * 100, 'LineWidth', 2);
end
yline(100, '-r', 'LineWidth', 1.4);
xlabel('Gait cycle (%)'); ylabel('% of speed limit'); title('Right joint speed at design speed'); grid on; xlim([0 100]);
legend(J, 'Location', 'best');
sgtitle('Joint trajectories and forward kinematics of the robot');
saveas(f2, [cfg.figPrefix '_joint_foot.png']);

f3 = figure('Name', 'Kinematic Analysis', 'Color', 'w', 'Position', [80 60 1500 450]);
cc = [0.2 0.6 0.2; 0 0 0; 0.8 0.2 0.2];
subplot(1,3,1); yyaxis left; hold on;
for s = 1:2
    plot(x, res(s).detJ, 'Color', colS{s}, 'LineWidth', 2);
end
ylabel('det J (mm^2)');
yyaxis right;
for s = 1:2
    plot(x, res(s).condJ, '--', 'Color', colS{s}, 'LineWidth', 1.3);
end
ylabel('cond(J)'); set(gca, 'YScale', 'log');
xlabel('Gait cycle (%)'); title('Singularity: det J (solid), cond J (dashed)'); grid on; xlim([0 100]);
subplot(1,3,2); hold on;
for c = 1:3
    plot(x, 100*adj(c).reach, 'Color', cc(c,:), 'LineWidth', 2);
end
yline(100, '-r', 'LineWidth', 1.4);
xlabel('Gait cycle (%)'); ylabel('% of leg reach'); title('Reach needed for the human path (right)'); grid on; xlim([0 100]);
legend(lensName, 'Location', 'best');
subplot(1,3,3); hold on;
for c = 1:3
    plot(adj(c).Ankle(:,1), adj(c).Ankle(:,2), 'Color', cc(c,:), 'LineWidth', 2);
end
axis equal; grid on; xlabel('x (mm)'); ylabel('z (mm)');
title('Ankle path, same joint angles, different link length');
legend(lensName, 'Location', 'best');
sgtitle('Kinematic analysis of the 3-DOF wearable robot');
saveas(f3, [cfg.figPrefix '_analysis.png']);

% ---------------- KINEMATIC EXPORT ----------------
if isfile(cfg.outFile), delete(cfg.outFile); end
for s = 1:2
    m = res(s).m;
    T = table(x, res(s).t, res(s).q(:,1), res(s).q(:,2), res(s).q(:,3), ...
        m.Hip(:,1), m.Hip(:,2), m.Knee(:,1), m.Knee(:,2), m.Ankle(:,1), m.Ankle(:,2), ...
        m.Heel(:,1), m.Heel(:,2), m.Toe(:,1), m.Toe(:,2), rad2deg(m.pitch), res(s).spJ, res(s).detJ, res(s).condJ, ...
        'VariableNames', {'GaitCycle_percent','Time_s','qHip_deg','qKnee_deg','qAnkle_deg', ...
        'HipX_mm','HipZ_mm','KneeX_mm','KneeZ_mm','AnkleX_mm','AnkleZ_mm','HeelX_mm','HeelZ_mm','ToeX_mm','ToeZ_mm', ...
        'FootAngle_deg','AnkleSpeed_mm_s','DetJ_mm2','CondJ'});
    writetable(T, cfg.outFile, 'Sheet', ['FK_' side{s}]);
end
sumRows = {
    'ThighLength', L1, 'mm'; 'ShankLength', L2, 'mm'; 'AnkleHeight', hA, 'mm'; 'FootPlateLength', Lf, 'mm'; 'HipHeight', hipZ, 'mm';
    'MaxReach_Left', 100*max(res(1).reach), '%'; 'MaxReach_Right', 100*max(res(2).reach), '%';
    'MinKneeFlexion_Left', min(res(1).q(:,2)), 'deg'; 'MinKneeFlexion_Right', min(res(2).q(:,2)), 'deg';
    'MaxFKIKError', max(cellfun(@(e) max(e(:)), {res.ikErr})), 'deg';
    'SoftMargin_Hip_min', min(softMargin(1,:)), 'deg'; 'SoftMargin_Knee_min', min(softMargin(2,:)), 'deg'; 'SoftMargin_Ankle_min', min(softMargin(3,:)), 'deg';
    'HardMargin_Hip_min', min(hardMargin(1,:)), 'deg'; 'HardMargin_Knee_min', min(hardMargin(2,:)), 'deg'; 'HardMargin_Ankle_min', min(hardMargin(3,:)), 'deg';
    'SpeedUsage_Hip_max', 100*max(velUse(1,:)), '%'; 'SpeedUsage_Knee_max', 100*max(velUse(2,:)), '%'; 'SpeedUsage_Ankle_max', 100*max(velUse(3,:)), '%';
    'AnkleShift_MinLinks', max(sqrt(sum((adj(1).Ankle - mNom.Ankle).^2, 2))), 'mm';
    'AnkleShift_MaxLinks', max(sqrt(sum((adj(3).Ankle - mNom.Ankle).^2, 2))), 'mm';
    'VerificationPassed', double(allOK), '1 = yes'};
writetable(cell2table(sumRows, 'VariableNames', {'Parameter','Value','Unit'}), cfg.outFile, 'Sheet', 'Summary');
fprintf('\nSaved: %s, %s_leg_movement.png, %s_joint_foot.png, %s_analysis.png\n', cfg.outFile, cfg.figPrefix, cfg.figPrefix, cfg.figPrefix);


%% ======================================================================
%%  PHASE III - ROBOT DYNAMIC MODEL
%% ======================================================================
dc.includeWearer = true;
dc.useGRF        = true;
dc.stancePct     = 60;
dc.g             = 9.81;
dc.bodyMass      = getPar(Dsn, 'WearerMass_kg', 70);
dc.nHarm         = 15;
dc.energyTol     = 0.01;
dc.gravTol       = 1e-6;
dc.margin        = 0.20;
dc.figPrefix     = 'robot_dynamics';

rmsf = @(v) sqrt(mean(v.^2, 1));

% ---------------- STEP 12 : MASS PROPERTIES ----------------
Lm  = [L1 L2 Lf] / 1000;
rm  = [getPar(Dsn,'ThighMass_kg',2.5) getPar(Dsn,'ShankMass_kg',1.8) getPar(Dsn,'FootMass_kg',0.9)];
rcr = [getPar(Dsn,'ThighCOM_ratio',0.45) getPar(Dsn,'ShankCOM_ratio',0.45)];
fx  = getPar(Dsn, 'FootCOM_x_mm', (0.5 - af)*Lf) / 1000;
fz  = getPar(Dsn, 'FootCOM_z_mm', -hA/2) / 1000;
hm  = dc.includeWearer * dc.bodyMass * [0.100 0.0465 0.0145];

[P,  mT, IT, cT] = buildP(rm, rcr, hm,       Lm, fx, fz, af, hA, dc.g, hipZ);
[P0, ~,  ~,  ~ ] = buildP(rm, rcr, [0 0 0],  Lm, fx, fz, af, hA, dc.g, hipZ);

linkName = {'Thigh','Shank','Foot'};
fprintf('\n========================================\n');
fprintf('   PHASE III - ROBOT MASS PROPERTIES\n');
fprintf('========================================\n');
fprintf('Wearer limb included : %s | wearer mass %.1f kg\n', yesNo(dc.includeWearer), dc.bodyMass);
fprintf('%-7s %10s %10s %10s %10s %18s %16s\n', 'Link', 'Length mm', 'Robot kg', 'Wearer kg', 'Total kg', 'COM from joint mm', 'I about COM kgm2');
for i = 1:3
    if i < 3
        comStr = sprintf('%.1f', 1000*cT(i));
    else
        comStr = sprintf('(%.1f, %.1f)', 1000*fx, 1000*fz);
    end
    fprintf('%-7s %10.1f %10.3f %10.3f %10.3f %18s %16.5f\n', linkName{i}, 1000*Lm(i), rm(i), hm(i), mT(i), comStr, IT(i));
end
fprintf('Total moving mass     : %.3f kg (robot only %.3f kg)\n', sum(mT), sum(rm));
fprintf('Foot COM is (x along foot, z normal to foot) relative to the ankle.\n');

% ---------------- STEP 13 : DYNAMIC MODEL ----------------
fprintf('\n========================================\n');
fprintf('   PHASE III - DYNAMIC MODEL\n');
fprintf('========================================\n');
fprintf('tau = M(q)*qdd + C(q,qd)*qd + G(q) + tau_ext   (tau = [hip knee ankle], Nm)\n');
fprintf('  M(q)   : sum_i [ m_i*Jv_i''*Jv_i + I_i*Jw_i''*Jw_i ]\n');
fprintf('  C(q,qd): Christoffel symbols of M(q)\n');
fprintf('  G(q)   : d/dq of sum_i m_i*g*z_i\n');
fprintf('  tau_ext: -Jc''*F_ground (vertical GRF in stance, centre of pressure heel to toe)\n');
fprintf('Joint angles are fitted with a %d-harmonic periodic Fourier series, so q, qd and qdd are mutually consistent.\n', dc.nHarm);
fprintf('Cases: FULL = robot + wearer + GRF | ROBOT = robot links only | HUMAN = FULL - ROBOT | ASSIST(a) = ROBOT + a*HUMAN\n');

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
hS   = 1e-6;
dPhi = 1e-5;
omDes = 2*pi/designCycle;
phi   = 2*pi*x/100;
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

    minEig = inf; symErr = 0; gravErr = 0;
    for i = 1:nP
        [Mi, ~, Gi] = dynTerms(qr(i,:), dqd(i,:), P);
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
    dyn(s).minEig = minEig; dyn(s).symErr = symErr; dyn(s).gravErr = gravErr; dyn(s).enErr = enErr;
    dyn(s).ddq = rad2deg(ddqd); dyn(s).dq = rad2deg(dqd); dyn(s).qf = qf;
    dyn(s).fitMax = fitMax; dyn(s).fitRms = fitRms; dyn(s).velDiff = velDiff;
end

fprintf('\n========================================\n');
fprintf('   JOINT-ANGLE FIT QUALITY (Fourier, %d harmonics)\n', dc.nHarm);
fprintf('========================================\n');
fprintf('%-6s %-6s %14s %14s %22s\n', 'Joint', 'Side', 'Max err (deg)', 'RMS err (deg)', 'dq vs recorded dq (%)');
for s = 1:2
    for j = 1:3
        fprintf('%-6s %-6s %14.3f %14.3f %22.1f\n', J{j}, side{s}, dyn(s).fitMax(j), dyn(s).fitRms(j), dyn(s).velDiff(j));
    end
end

fprintf('\n========================================\n');
fprintf('   CASE FULL: REQUIRED JOINT TORQUE (robot + wearer + GRF, design speed)\n');
fprintf('========================================\n');
fprintf('%-6s %-6s %10s %10s %11s %10s %12s %11s %9s\n', 'Joint', 'Side', 'Peak Nm', 'RMS Nm', 'Mean|t| Nm', 'Peak rec.', 'Peak P (W)', 'Peak limit', 'Usage');
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
        if isfinite(pkLim(j))
            limStr = sprintf('%.1f Nm', pkLim(j)); useStr = sprintf('%.0f %%', 100*pkT(j,s)/pkLim(j));
        else
            limStr = 'n/a'; useStr = 'n/a';
        end
        fprintf('%-6s %-6s %10.2f %10.2f %11.2f %10.2f %12.1f %11s %9s\n', J{j}, side{s}, pkT(j,s), rmsT(j,s), mean(abs(tj)), pkRec(j,s), pkP(j,s), limStr, useStr);
    end
end
fprintf('Peak limit is read from a torque column of the Requirements/Actuators sheets.\n');

fprintf('\n========================================\n');
fprintf('   CASE ROBOT: TORQUE TO MOVE THE ROBOT LINKS ONLY (transparent mode)\n');
fprintf('========================================\n');
fprintf('%-6s %-6s %12s %12s %11s %9s\n', 'Joint', 'Side', 'Peak Nm', 'RMS Nm', 'Peak limit', 'Usage');
for s = 1:2
    for j = 1:3
        if isfinite(pkLim(j))
            limStr = sprintf('%.1f Nm', pkLim(j)); useStr = sprintf('%.0f %%', 100*pkR(j,s)/pkLim(j));
        else
            limStr = 'n/a'; useStr = 'n/a';
        end
        fprintf('%-6s %-6s %12.2f %12.2f %11s %9s\n', J{j}, side{s}, pkR(j,s), rmsf(dyn(s).tauR(:,j)), limStr, useStr);
    end
end

fprintf('\n========================================\n');
fprintf('   TORQUE COMPONENTS, FULL CASE (RMS, Nm)\n');
fprintf('========================================\n');
fprintf('%-6s %-6s %12s %12s %12s %12s %12s\n', 'Joint', 'Side', 'M*qdd', 'C*qd', 'G', 'tau_ext', 'Total');
for s = 1:2
    for j = 1:3
        fprintf('%-6s %-6s %12.2f %12.2f %12.2f %12.2f %12.2f\n', J{j}, side{s}, rmsf(dyn(s).Tin(:,j)), rmsf(dyn(s).Tcor(:,j)), ...
                rmsf(dyn(s).Tg(:,j)), rmsf(dyn(s).Text(:,j)), rmsT(j,s));
    end
end

% ---------------- ACTUATOR CAPACITY / ASSIST ANALYSIS ----------------
alG = linspace(0, 1, 2001);
alphaMax = nan(3,2);
for s = 1:2
    for j = 1:3
        if ~isfinite(pkLim(j)), continue; end
        okv = false(size(alG));
        for a = 1:numel(alG)
            okv(a) = max(abs(dyn(s).tauR(:,j) + alG(a)*dyn(s).tauH(:,j))) <= pkLim(j);
        end
        if any(okv), alphaMax(j,s) = alG(find(okv, 1, 'last')); end
    end
end
alphaJ = min(alphaMax, [], 2)';
reqRating = (1 + dc.margin) * max(pkT, [], 2)';

fprintf('\n========================================\n');
fprintf('   ACTUATOR CAPACITY AND ASSIST LEVEL (worst of both legs)\n');
fprintf('========================================\n');
fprintf('%-6s %11s %14s %14s %16s %18s %16s\n', 'Joint', 'Limit Nm', 'Full peak Nm', 'Robot peak Nm', 'Full support', 'Max assist ratio', 'Rating needed*');
for j = 1:3
    fp = max(pkT(j,:)); rp = max(pkR(j,:));
    if isfinite(pkLim(j))
        if fp <= pkLim(j), fs = 'OK'; else, fs = sprintf('SHORT %.0f Nm', fp - pkLim(j)); end
        if isnan(alphaJ(j)), as = 'none'; else, as = sprintf('%.0f %%', 100*alphaJ(j)); end
        fprintf('%-6s %11.1f %14.2f %14.2f %16s %18s %13.1f Nm\n', J{j}, pkLim(j), fp, rp, fs, as, reqRating(j));
    else
        fprintf('%-6s %11s %14.2f %14.2f %16s %18s %13.1f Nm\n', J{j}, 'n/a', fp, rp, 'n/a', 'n/a', reqRating(j));
    end
end
fprintf('* peak of the FULL case plus %.0f %% design margin.\n', 100*dc.margin);
fprintf('Max assist ratio = largest share of the human joint torque (wearer limb + ground load) the actuator can supply,\n');
fprintf('on top of moving its own links, without exceeding the peak-torque limit.\n');

fprintf('\n========================================\n');
fprintf('   JOINT POWER AND WORK PER CYCLE, FULL CASE (design speed)\n');
fprintf('========================================\n');
fprintf('%-6s %-6s %12s %12s %12s %12s\n', 'Joint', 'Side', 'Pos work J', 'Neg work J', 'Net work J', 'Avg |P| W');
for s = 1:2
    for j = 1:3
        fprintf('%-6s %-6s %12.2f %12.2f %12.2f %12.2f\n', J{j}, side{s}, dyn(s).Wpos(j), dyn(s).Wneg(j), dyn(s).Wpos(j)+dyn(s).Wneg(j), mean(abs(dyn(s).pw(:,j))));
    end
end

fprintf('\n========================================\n');
fprintf('   FULL-CASE TORQUE AT KEY GAIT EVENTS (right leg, Nm)\n');
fprintf('========================================\n');
evP = [0 15 30 50 60 75 100];
fprintf('%-12s %10s %10s %10s\n', 'Gait cycle', 'Hip', 'Knee', 'Ankle');
for e = 1:numel(evP)
    [~, k] = min(abs(x - evP(e)));
    fprintf('%-10.0f %% %10.2f %10.2f %10.2f\n', x(k), dyn(2).tau(k,1), dyn(2).tau(k,2), dyn(2).tau(k,3));
end

symOK  = max([dyn.symErr]) < 1e-9;
pdOK   = min([dyn.minEig]) > 0;
gravOK = max([dyn.gravErr]) < dc.gravTol;
enOK   = max([dyn.enErr]) < dc.energyTol;
dynOK  = symOK && pdOK && gravOK && enOK && all(isfinite(dyn(1).tau(:))) && all(isfinite(dyn(2).tau(:)));
limKnown = isfinite(pkLim);
fullOK  = ~any(limKnown) || all(all(pkT(limKnown,:)  <= repmat(pkLim(limKnown)', 1, 2)));
robotOK = ~any(limKnown) || all(all(pkR(limKnown,:)  <= repmat(pkLim(limKnown)', 1, 2)));

fprintf('\n========================================\n');
fprintf('   DYNAMIC MODEL VERIFICATION\n');
fprintf('========================================\n');
fprintf('M(q) symmetric (max asym %.1e)               : %s\n', max([dyn.symErr]), passStr(symOK));
fprintf('M(q) positive definite (min eig %.4f)        : %s\n', min([dyn.minEig]), passStr(pdOK));
fprintf('G(q) vs numerical dV/dq (max err %.1e Nm)   : %s\n', max([dyn.gravErr]), passStr(gravOK));
fprintf('Power balance qd''*tau = dE/dt (max err %.3f%%) : %s\n', 100*max([dyn.enErr]), passStr(enOK));
if dynOK
    fprintf('The dynamic model is consistent.\n');
else
    fprintf('At least one model check failed. Review the tables above.\n');
end

fprintf('\n========================================\n');
fprintf('   ACTUATOR SUFFICIENCY\n');
fprintf('========================================\n');
if any(limKnown)
    fprintf('Full support of the gait (robot + wearer + GRF) : %s\n', passStr(fullOK));
    fprintf('Moving the robot links only                     : %s\n', passStr(robotOK));
    for j = 1:3
        if isnan(alphaJ(j))
            fprintf('  %-6s : actuator cannot supply the required torque even in transparent mode\n', J{j});
        else
            fprintf('  %-6s : actuator covers up to %.0f %% of the human joint torque\n', J{j}, 100*alphaJ(j));
        end
    end
    if ~fullOK
        fprintf('The current actuators cannot fully replace the joint torque in the %.2f s design cycle.\n', designCycle);
        fprintf('Use them in assist mode at the ratios above, or raise the ratings to the "Rating needed" column.\n');
    end
else
    fprintf('No torque rating found in the input sheets, so actuator sufficiency was not evaluated.\n');
end

% ---------------- PHASE III PLOTS ----------------
f5 = figure('Name', 'Robot Joint Torque', 'Color', 'w', 'Position', [40 30 1500 780]);
for j = 1:3
    subplot(2,3,j); hold on;
    for s = 1:2
        plot(x, dyn(s).tau(:,j), 'Color', colS{s}, 'LineWidth', 2);
    end
    if isfinite(pkLim(j)), yline(pkLim(j), '-r', 'LineWidth', 1.3); yline(-pkLim(j), '-r', 'LineWidth', 1.3); end
    xlabel('Gait cycle (%)'); ylabel('Nm'); title([J{j} ' torque, full case (design speed)']); grid on; xlim([0 100]);
    if j == 1, legend('Left','Right','Location','best'); end
end
cmpCol = [0 0.447 0.741; 0.929 0.694 0.125; 0.466 0.674 0.188; 0.494 0.184 0.556];
for j = 1:3
    subplot(2,3,3+j); hold on;
    plot(x, dyn(2).Tin(:,j),  'Color', cmpCol(1,:), 'LineWidth', 1.5);
    plot(x, dyn(2).Tcor(:,j), 'Color', cmpCol(2,:), 'LineWidth', 1.5);
    plot(x, dyn(2).Tg(:,j),   'Color', cmpCol(3,:), 'LineWidth', 1.5);
    plot(x, dyn(2).Text(:,j), 'Color', cmpCol(4,:), 'LineWidth', 1.5);
    plot(x, dyn(2).tau(:,j),  'k', 'LineWidth', 2.2);
    xlabel('Gait cycle (%)'); ylabel('Nm'); title([J{j} ' torque components (right)']); grid on; xlim([0 100]);
    if j == 1, legend('M*qdd','C*qd','G','tau ext','Total','Location','best'); end
end
sgtitle('Robot dynamic model: required joint torque and its components');
saveas(f5, [dc.figPrefix '_torque.png']);

f6 = figure('Name', 'Robot Dynamics Checks', 'Color', 'w', 'Position', [60 60 1300 450]);
subplot(1,2,1); yyaxis left;
plot(x, Fz, 'LineWidth', 2); ylabel('Vertical GRF (N)');
yyaxis right;
plot(x, 100*frac.*stc, '--', 'LineWidth', 1.5); ylabel('Centre of pressure (% foot from heel)');
xlabel('Gait cycle (%)'); title('External load model (tau ext)'); grid on; xlim([0 100]);
subplot(1,2,2); hold on;
plot(x, dyn(2).Pn, 'Color', colS{2}, 'LineWidth', 2);
plot(x, dyn(2).dE, 'k:', 'LineWidth', 1.8);
xlabel('Gait cycle (%)'); ylabel('W'); grid on; xlim([0 100]);
title(sprintf('Energy check (right): max error %.3f %%', 100*dyn(2).enErr));
legend('qd''*(M qdd + C qd + G)','dE/dt','Location','best');
sgtitle('Dynamic model load and verification');
saveas(f6, [dc.figPrefix '_analysis.png']);

% ---------------- PHASE III EXPORT ----------------
dynNames = {'GaitCycle_percent','Time_s','TauHip_Nm','TauKnee_Nm','TauAnkle_Nm', ...
            'InertialHip_Nm','InertialKnee_Nm','InertialAnkle_Nm', ...
            'CoriolisHip_Nm','CoriolisKnee_Nm','CoriolisAnkle_Nm', ...
            'GravityHip_Nm','GravityKnee_Nm','GravityAnkle_Nm', ...
            'ExtHip_Nm','ExtKnee_Nm','ExtAnkle_Nm', ...
            'PowerHip_W','PowerKnee_W','PowerAnkle_W','AccHip_deg_s2','AccKnee_deg_s2','AccAnkle_deg_s2','GRF_N', ...
            'RobotOnlyHip_Nm','RobotOnlyKnee_Nm','RobotOnlyAnkle_Nm','HumanShareHip_Nm','HumanShareKnee_Nm','HumanShareAnkle_Nm'};
for s = 1:2
    D = [x, td, dyn(s).tau, dyn(s).Tin, dyn(s).Tcor, dyn(s).Tg, dyn(s).Text, dyn(s).pw, dyn(s).ddq, Fz, dyn(s).tauR, dyn(s).tauH];
    writetable(array2table(D, 'VariableNames', dynNames), cfg.outFile, 'Sheet', ['Dyn_' side{s}]);
end
massTab = table(linkName', 1000*Lm', rm', hm', mT', [1000*cT 1000*fx]', IT', ...
    'VariableNames', {'Link','Length_mm','RobotMass_kg','WearerMass_kg','TotalMass_kg','COM_mm','Inertia_kgm2'});
writetable(massTab, cfg.outFile, 'Sheet', 'MassProps');
assistTab = table(J', pkLim', max(pkT,[],2), max(pkR,[],2), alphaJ', reqRating', ...
    'VariableNames', {'Joint','PeakLimit_Nm','FullPeak_Nm','RobotOnlyPeak_Nm','MaxAssistRatio','RequiredRating_Nm'});
writetable(assistTab, cfg.outFile, 'Sheet', 'ActuatorCapacity');
dynRows = {};
for s = 1:2
    for j = 1:3
        dynRows(end+1,:) = {sprintf('PeakTorque_%s_%s', J{j}, side{s}), pkT(j,s), 'Nm'};
        dynRows(end+1,:) = {sprintf('RMSTorque_%s_%s',  J{j}, side{s}), rmsT(j,s), 'Nm'};
        dynRows(end+1,:) = {sprintf('PeakPower_%s_%s',  J{j}, side{s}), pkP(j,s), 'W'};
    end
end
dynRows(end+1,:) = {'PowerBalanceError_max', 100*max([dyn.enErr]), '%'};
dynRows(end+1,:) = {'GravityVectorError_max', max([dyn.gravErr]), 'Nm'};
dynRows(end+1,:) = {'DynamicsVerificationPassed', double(dynOK), '1 = yes'};
dynRows(end+1,:) = {'FullSupportWithinActuators', double(fullOK), '1 = yes'};
writetable(cell2table(dynRows, 'VariableNames', {'Parameter','Value','Unit'}), cfg.outFile, 'Sheet', 'DynSummary');
fprintf('\nSaved: %s (Dyn_Left, Dyn_Right, MassProps, ActuatorCapacity, DynSummary), %s_torque.png, %s_analysis.png\n', ...
        cfg.outFile, dc.figPrefix, dc.figPrefix);


%% ======================================================================
%%  STEP 14 : ACTUATOR REQUIREMENTS AND SELECTION
%% ======================================================================
ac.SF_peak     = 1.5;
ac.SF_cont     = 1.3;
ac.SF_speed    = 1.2;
ac.SF_power    = 1.3;
ac.assistRatio = 0.5;
ac.designCase  = 1;       % 1 FULL, 2 ASSIST, 3 ROBOT
ac.gear        = [getPar(Dsn,'HipGearRatio',100) getPar(Dsn,'KneeGearRatio',100) getPar(Dsn,'AnkleGearRatio',100)];
ac.eta         = getPar(Dsn, 'GearEfficiency', 0.85);
ac.figPrefix   = 'robot_actuator';

caseName = {'FULL','ASSIST','ROBOT'};
caseDesc = {'FULL support (robot + wearer + ground load)', ...
            sprintf('ASSIST %.0f %% of human joint torque + robot links', 100*ac.assistRatio), ...
            'ROBOT links only (transparent mode)'};

fprintf('\n========================================\n');
fprintf('   STEP 14 - ACTUATOR REQUIREMENTS\n');
fprintf('========================================\n');
fprintf('Source        : dynamic simulation at the %.2f s design cycle (both legs, worst leg is used for sizing)\n', designCycle);
fprintf('Continuous T  : RMS joint torque over one gait cycle (thermal sizing basis)\n');
fprintf('Joint power   : p = tau * omega (W); omega in rad/s\n');
fprintf('Selection     : T_select = SF x T_required, rounded up to the next R10 standard size\n');
fprintf('Comparison    : required torque / existing actuator limit (>100 %% = limit exceeded)\n');
fprintf('Motor torque  : standard size / (gear ratio x efficiency)\n');
fprintf('Safety factors: peak torque x%.2f | continuous torque x%.2f | speed x%.2f | power x%.2f\n', ac.SF_peak, ac.SF_cont, ac.SF_speed, ac.SF_power);
fprintf('Gear train    : N = [%g %g %g], efficiency %.2f\n', ac.gear(1), ac.gear(2), ac.gear(3), ac.eta);

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
Tcont = Trms;

TselPk = ac.SF_peak  * Tpk;
TselC  = ac.SF_cont  * Tcont;
WselR  = ac.SF_speed * Wpk;
Psel   = ac.SF_power * Ppk;
stdPk = zeros(3,nC); stdC = zeros(3,nC); motT = zeros(3,nC); motTc = zeros(3,nC); motRpm = zeros(3,nC);
for c = 1:nC
    for j = 1:3
        stdPk(j,c)  = stdSize(TselPk(j,c));
        stdC(j,c)   = stdSize(TselC(j,c));
        motT(j,c)   = stdPk(j,c) / (ac.gear(j) * ac.eta);
        motTc(j,c)  = stdC(j,c)  / (ac.gear(j) * ac.eta);
        motRpm(j,c) = WselR(j,c) * ac.gear(j) * 60 / (2*pi);
    end
end

for c = 1:nC
    fprintf('\n----------------------------------------\n');
    fprintf(' CASE %s : %s\n', caseName{c}, caseDesc{c});
    fprintf('----------------------------------------\n');
    fprintf('REQUIRED (from simulation)\n');
    fprintf('%-6s %9s %9s %11s %10s %10s %10s %10s %9s %10s\n', 'Joint', 'Peak T', 'RMS T', 'Contin. T', 'Peak w', 'Peak w', 'RMS w', 'Peak w', 'Peak P', 'RMS P');
    fprintf('%-6s %9s %9s %11s %10s %10s %10s %10s %9s %10s\n', '', 'Nm', 'Nm', 'Nm', 'deg/s', 'rad/s', 'deg/s', 'rpm', 'W', 'W');
    for j = 1:3
        fprintf('%-6s %9.2f %9.2f %11.2f %10.1f %10.2f %10.1f %10.1f %9.1f %10.1f\n', J{j}, Tpk(j,c), Trms(j,c), Tcont(j,c), ...
                rad2deg(Wpk(j,c)), Wpk(j,c), rad2deg(Wrms(j,c)), Wpk(j,c)*60/(2*pi), Ppk(j,c), Prms(j,c));
    end
    fprintf('Average positive (driving) power per cycle, W:  Hip %.1f | Knee %.1f | Ankle %.1f\n', Ppos(1,c), Ppos(2,c), Ppos(3,c));
    fprintf('SELECTION (T_select = SF x T_required)\n');
    fprintf('%-6s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n', 'Joint', 'Peak T sel', 'Rated T sel', 'Speed sel', 'Power sel', 'Std peak', 'Std rated', 'Motor T pk', 'Motor T rtd', 'Motor speed');
    fprintf('%-6s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n', '', 'Nm', 'Nm', 'deg/s', 'W', 'Nm', 'Nm', 'Nm', 'Nm', 'rpm');
    for j = 1:3
        fprintf('%-6s %12.2f %12.2f %12.1f %12.1f %12.2f %12.2f %12.3f %12.3f %12.0f\n', J{j}, TselPk(j,c), TselC(j,c), rad2deg(WselR(j,c)), Psel(j,c), ...
                stdPk(j,c), stdC(j,c), motT(j,c), motTc(j,c), motRpm(j,c));
    end
    fprintf('COMPARISON WITH THE ACTUATORS IN THE DESIGN SHEET (required / limit)\n');
    fprintf('%-6s %14s %14s %14s %14s %14s %14s\n', 'Joint', 'Peak limit Nm', 'Peak use', 'Rated limit Nm', 'Rated use', 'Speed limit', 'Speed status');
    for j = 1:3
        fprintf('%-6s %14s %14s %14s %14s %14s %14s\n', J{j}, numStr(pkLim(j), '%.1f'), utilStr(pkLim(j), Tpk(j,c)), ...
                numStr(rtLim(j), '%.1f'), utilStr(rtLim(j), Tcont(j,c)), numStr(velLim(j), '%.0f d/s'), statStr(velLim(j), rad2deg(WselR(j,c))));
    end
end

dcs = ac.designCase;
fprintf('\n========================================\n');
fprintf('   ACTUATOR SPECIFICATION TO PROCURE (case %s)\n', caseName{dcs});
fprintf('========================================\n');
fprintf('%-6s %16s %16s %16s %14s %16s %16s %16s\n', 'Joint', 'Peak torque Nm', 'Rated torque Nm', 'No-load speed', 'Power W', 'Gear ratio N', 'Motor peak Nm', 'Motor rated Nm');
fprintf('%-6s %16s %16s %16s %14s %16s %16s %16s\n', '', '(>=)', '(>=)', 'deg/s (>=)', '(>=)', '', '(>=)', '(>=)');
for j = 1:3
    fprintf('%-6s %16.1f %16.1f %16.0f %14.0f %16g %16.3f %16.3f\n', J{j}, stdPk(j,dcs), stdC(j,dcs), ceil(rad2deg(WselR(j,dcs))), ceil(Psel(j,dcs)), ac.gear(j), motT(j,dcs), motTc(j,dcs));
end
fprintf('Peak torque = %.2f x peak required, rated torque = %.2f x RMS required, both rounded up to the next R10 standard value.\n', ac.SF_peak, ac.SF_cont);
fprintf('Peak torque must be available at the joint for at least the stance burst; rated torque must be sustained continuously.\n');

fprintf('\n========================================\n');
fprintf('   STEP 14 CHECKS\n');
fprintf('========================================\n');
chkSF  = all([ac.SF_peak ac.SF_cont ac.SF_speed ac.SF_power] >= 1);
chkSel = all(TselPk(:) >= Tpk(:)) && all(TselC(:) >= Tcont(:)) && all(WselR(:) >= Wpk(:)) && all(Psel(:) >= Ppk(:));
chkStd = all(stdPk(:) >= TselPk(:) - 1e-9) && all(stdC(:) >= TselC(:) - 1e-9);
chkPk  = all(Tpk(:) >= Trms(:) - 1e-9);
chkMat = max(abs(Tpk(:,1) - max(pkT,[],2))) < 1e-9;
chkMot = all(abs(motT(:)  - stdPk(:) ./ repmat(ac.gear(:)*ac.eta, nC, 1)) < 1e-12) && ...
         all(abs(motTc(:) - stdC(:)  ./ repmat(ac.gear(:)*ac.eta, nC, 1)) < 1e-12);
fprintf('Safety factors >= 1                              : %s\n', passStr(chkSF));
fprintf('Selected values >= required values               : %s\n', passStr(chkSel));
fprintf('Standard sizes >= selected values                : %s\n', passStr(chkStd));
fprintf('Peak torque >= RMS torque                        : %s\n', passStr(chkPk));
fprintf('Step 14 FULL peak torque equals Phase III result : %s\n', passStr(chkMat));
fprintf('Motor torque = standard size / (N x eta)         : %s\n', passStr(chkMot));

% ---------------- STEP 14 PLOTS ----------------
f8 = figure('Name', 'Actuator Requirements', 'Color', 'w', 'Position', [30 20 1500 780]);
for j = 1:3
    subplot(2,3,j); hold on;
    hh = gobjects(1,4);
    for s = 1:2
        hh(s) = plot(dyn(s).dq(:,j), tauCase{1,s}(:,j), '.-', 'Color', colS{s}, 'LineWidth', 1.2, 'MarkerSize', 8);
    end
    for s = 1:2
        hh(2+s) = plot(dyn(s).dq(:,j), tauCase{3,s}(:,j), '.', 'Color', colS{s}*0.5 + 0.5, 'MarkerSize', 6);
    end
    ws = rad2deg(WselR(j,dcs)); ts = TselPk(j,dcs);
    rectangle('Position', [-ws -ts 2*ws 2*ts], 'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1.6);
    if isfinite(pkLim(j)) && isfinite(velLim(j))
        rectangle('Position', [-velLim(j) -pkLim(j) 2*velLim(j) 2*pkLim(j)], 'EdgeColor', [0 0.6 0], 'LineWidth', 1.6);
    end
    yline(0, 'k-'); xline(0, 'k-');
    xlabel('Joint speed (deg/s)'); ylabel('Torque (Nm)'); grid on;
    title([J{j} ' torque-speed (red: selected, green: existing)']);
    if j == 1, safeLegend(hh, {'Left full','Right full','Left robot','Right robot'}, 'northeast'); end
end
for j = 1:3
    subplot(2,3,3+j); hold on;
    hp = gobjects(1,3);
    for s = 1:2
        hp(s) = plot(x, pwCase{1,s}(:,j), 'Color', colS{s}, 'LineWidth', 2);
    end
    hp(3) = yline(Psel(j,dcs), '--r', 'LineWidth', 1.3); yline(-Psel(j,dcs), '--r', 'LineWidth', 1.3);
    yline(0, 'k-');
    xlabel('Gait cycle (%)'); ylabel('W'); grid on; xlim([0 100]);
    title(sprintf('%s power, selected %.0f W', J{j}, Psel(j,dcs)));
    if j == 1, safeLegend(hp, {'Left','Right','Selected power'}, 'northeast'); end
end
sgtitle(sprintf('Step 14 actuator requirements (design cycle %.2f s, selection case %s)', designCycle, caseName{dcs}));
saveas(f8, [ac.figPrefix '_requirements.png']);

% ---------------- STEP 14 EXPORT ----------------
actRows = {};
for c = 1:nC
    for j = 1:3
        actRows(end+1,:) = {caseName{c}, J{j}, Tpk(j,c), Trms(j,c), Tcont(j,c), rad2deg(Wpk(j,c)), rad2deg(Wrms(j,c)), ...
                            Ppk(j,c), Prms(j,c), Ppos(j,c), TselPk(j,c), TselC(j,c), rad2deg(WselR(j,c)), Psel(j,c), ...
                            stdPk(j,c), stdC(j,c), ac.gear(j), motT(j,c), motTc(j,c), motRpm(j,c), pkLim(j), rtLim(j), velLim(j), ...
                            100*Tpk(j,c)/pkLim(j), 100*Tcont(j,c)/rtLim(j)};
    end
end
actTab = cell2table(actRows, 'VariableNames', {'Case','Joint','PeakTorque_Nm','RMSTorque_Nm','ContinuousTorque_Nm', ...
    'PeakSpeed_deg_s','RMSSpeed_deg_s','PeakPower_W','RMSPower_W','AvgPositivePower_W', ...
    'SelPeakTorque_Nm','SelRatedTorque_Nm','SelSpeed_deg_s','SelPower_W','StdPeakTorque_Nm','StdRatedTorque_Nm', ...
    'GearRatio','MotorPeakTorque_Nm','MotorRatedTorque_Nm','MotorSpeed_rpm','ExistingPeakLimit_Nm','ExistingRatedLimit_Nm','ExistingSpeedLimit_deg_s', ...
    'PeakUse_percent','RatedUse_percent'});
writetable(actTab, cfg.outFile, 'Sheet', 'ActuatorSelection');
sfRows = {'SF_PeakTorque', ac.SF_peak, '-'; 'SF_ContinuousTorque', ac.SF_cont, '-'; 'SF_Speed', ac.SF_speed, '-'; 'SF_Power', ac.SF_power, '-'; ...
          'AssistRatio', ac.assistRatio, '-'; 'GearEfficiency', ac.eta, '-'; 'DesignCycleTime', designCycle, 's'; ...
          'SelectionCase', double(dcs), '1 FULL, 2 ASSIST, 3 ROBOT'};
writetable(cell2table(sfRows, 'VariableNames', {'Parameter','Value','Unit'}), cfg.outFile, 'Sheet', 'SelectionInputs');
fprintf('\nSaved: %s (ActuatorSelection, SelectionInputs), %s_requirements.png\n', cfg.outFile, ac.figPrefix);

%% ============================= FUNCTIONS =============================
function safeLegend(h, labels, loc)
    try
        legend(h, labels, 'Location', loc);
    catch
        try
            legend(h, labels);
        catch
        end
    end
end

function v = getPar(T, name, def)
    v = T.Value(strcmp(string(T.Parameter), name));
    if isempty(v)
        if nargin < 3
            error('Parameter %s not found.', name);
        end
        v = def;
    end
    v = v(1);
end

function s = passStr(ok)
    if ok, s = 'PASS'; else, s = 'CHECK'; end
end

function s = yesNo(ok)
    if ok, s = 'yes'; else, s = 'no'; end
end

function s = numStr(v, fmt)
    if isfinite(v), s = sprintf(fmt, v); else, s = 'n/a'; end
end

function s = statStr(lim, need)
    if ~isfinite(lim)
        s = 'n/a';
    elseif lim >= need
        s = sprintf('OK (%.0f%%)', 100*need/lim);
    else
        s = sprintf('SHORT %.0f%%', 100*(need/lim - 1));
    end
end

function s = utilStr(lim, need)
    if ~isfinite(lim)
        s = 'n/a';
    elseif need <= lim
        s = sprintf('%.0f%% OK', 100*need/lim);
    else
        s = sprintf('%.0f%% OVER', 100*need/lim);
    end
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

function [P, mT, IT, cT] = buildP(rm, rcr, hm, Lm, fx, fz, af, hA, g, hipZ)
    rk  = 1/sqrt(12);
    hcr = [0.433 0.433];
    hrg = [0.323 0.302 0.475];
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

function [M, G] = massMat(q, P)
    Tm  = [1 0 0; 1 -1 0; 1 -1 1];
    th1 = q(1); th2 = q(1) - q(2); pit = q(3) + th2;
    d1  = [cos(th1); sin(th1)];
    d2  = [cos(th2); sin(th2)];
    u   = [cos(pit); sin(pit)];
    nv  = [-sin(pit); cos(pit)];
    Jv  = cell(1,3);
    Jv{1} = P.c1*d1*Tm(1,:);
    Jv{2} = P.L1*d1*Tm(1,:) + P.c2*d2*Tm(2,:);
    Jv{3} = P.L1*d1*Tm(1,:) + P.L2*d2*Tm(2,:) + (P.fx*nv - P.fz*u)*Tm(3,:);
    M = zeros(3); G = zeros(3,1);
    for i = 1:3
        M = M + P.m(i)*(Jv{i}'*Jv{i}) + P.I(i)*(Tm(i,:)'*Tm(i,:));
        G = G + P.g*P.m(i)*Jv{i}(2,:)';
    end
end

function [M, C, G] = dynTerms(q, dq, P)
    h  = 1e-6;
    [M, G] = massMat(q, P);
    dM = zeros(3,3,3);
    for k = 1:3
        e = zeros(1,3); e(k) = h;
        dM(:,:,k) = (massMat(q + e, P) - massMat(q - e, P)) / (2*h);
    end
    C = zeros(3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                C(i,j) = C(i,j) + 0.5*(dM(i,j,k) + dM(i,k,j) - dM(j,k,i))*dq(k);
            end
        end
    end
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