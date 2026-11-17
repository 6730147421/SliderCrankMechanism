function build_SliderCrank_Simscape(RPM, theta0_deg, nRev)
%BUILD_SLIDERCRANK_SIMSCAPE  Build + run the Simscape Multibody model of the
%   3เกลอ Crank–Slider mechanism (3D resin print, scale 300 %).
%
%   วิธีใช้ (Usage) — วางไฟล์นี้ไว้ในโฟลเดอร์เดียวกับโฟลเดอร์ Parts_300pct แล้วรัน:
%       build_SliderCrank_Simscape              % ค่าเริ่มต้น 30 RPM, theta0 = 90 deg, 2 รอบ
%       build_SliderCrank_Simscape(60, 0, 3)    % 60 RPM, เริ่มที่ 0 deg, 3 รอบ
%
%   สิ่งที่ได้:
%     1) ไฟล์ SliderCrank_3Kelor_Simscape.slx (Block diagram สำหรับหน้า report)
%     2) หน้าต่าง Mechanics Explorer แสดงการเคลื่อนที่ของชิ้นงาน 3D (STL จริง)
%     3) กราฟ Position / Velocity / Acceleration ของ slider เทียบ Simscape กับสมการ
%     4) SliderCrank_results.csv สำหรับเทียบกับ Mathematical / Physical model
%
%   ข้อมูลแบบ (scaled 300 %):  r = 45 mm, L = 90 mm, e = 0 mm, Stroke = 2r = 90 mm
%   Frame ~ 210 x 100.95 x 97.5 mm  — ระยะ link วัดแบบ center-to-center
%   (ตรวจสอบจากไฟล์ CAD: แกน crank (x,z) = (-30,15), Pin B = (-30,60),
%    Pin A = (47.94,15)  ->  r = 45.00, L = 90.00, e = 0.00 mm)
%
%   Requires: Simulink, Simscape, Simscape Multibody (R2019b or newer).

if nargin < 1 || isempty(RPM),        RPM = 30;        end
if nargin < 2 || isempty(theta0_deg), theta0_deg = 90; end   % 90 = ท่าเดียวกับไฟล์ CAD
if nargin < 3 || isempty(nRev),       nRev = 2;        end

%% ---------------- Parameters (mm, deg) ----------------
P.r   = 45;          % crank radius
P.L   = 90;          % connecting rod length (center-to-center)
P.e   = 0;           % offset (รอตรวจสอบ) — geometry STL ถูกวาดที่ e = 0
P.Ocad = [-30 15];   % crank axis in CAD coordinates (x, z) [mm]
P.rho = 1150;        % resin density [kg/m^3] (ใช้คำนวณมวลจาก geometry)

mdl  = 'SliderCrank_3Kelor_Simscape';
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
partsDir = findPartsDir(here);
useSTL   = ~isempty(partsDir);
if ~useSTL
    warning('SC:noSTL', ['ไม่พบโฟลเดอร์ Parts_300pct — จะใช้รูปทรงอย่างง่าย (Brick) แทน STL.\n' ...
        'วางไฟล์นี้ไว้ข้างโฟลเดอร์ Parts_300pct เพื่อให้เห็นชิ้นงานจริง']);
else
    fprintf('Using STL parts from: %s\n', partsDir);
end

%% ---------------- Geometry / frame bookkeeping ----------------
% CAD frame: x = along slider, y = depth, z = up.
% Model frame: X = CAD x, Y = CAD z (up), Z = -CAD y  -> mechanism in the XY plane,
% every revolute axis = model Z (Simscape joints rotate about Z).
Rc2m = [1 0 0; 0 0 1; 0 -1 0];
Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
Ry = @(a) [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];

O     = [P.Ocad(1); P.Ocad(2); 0];             % crank axis (model)
thCAD = 90;                                     % crank angle in the CAD pose
PB    = O + Rz(thCAD)*[P.r; 0; 0];              % crank pin  (-30, 60)
xAcad = P.r*cosd(thCAD) + sqrt(P.L^2 - (P.r*sind(thCAD))^2);
PA    = O + [xAcad; 0; 0];                      % slider pin (47.94, 15)
phCAD = atan2d(PA(2)-PB(2), PA(1)-PB(1));       % rod angle in CAD pose (-30 deg)

% Rigid transforms that place the STL (CAD coords) into each body frame
T.base.R   = Rc2m;                 T.base.t   = [0; 0; 0];
T.crank.R  = Rz(-thCAD)*Rc2m;      T.crank.t  = -Rz(-thCAD)*O;
T.rod.R    = Rz(-phCAD)*Rc2m;      T.rod.t    = -Rz(-phCAD)*PB;
T.slider.R = Rc2m;                 T.slider.t = -PA;

% Initial assembly guesses (low priority) so the solver picks the right branch
x0   = P.r*cosd(theta0_deg) + sqrt(P.L^2 - (P.r*sind(theta0_deg) - P.e)^2);
phi0 = atan2d(P.e - P.r*sind(theta0_deg), x0 - P.r*cosd(theta0_deg));
relB0 = wrapTo180local(phi0 - theta0_deg);      % crank -> rod joint angle

% Colours (same as the printed parts)
C.frame  = [0.70 0.72 0.73];  C.frame2 = [0.40 0.42 0.44];
C.white  = [0.97 0.97 0.97];  C.red    = [0.84 0.17 0.14];
C.pink   = [0.97 0.55 0.73];  C.pin    = [0.82 0.83 0.84];
C.shaft  = [0.79 0.80 0.81];

%% ---------------- New model ----------------
if bdIsLoaded(mdl), close_system(mdl, 0); end
slx = fullfile(here, [mdl '.slx']);
if exist(slx, 'file'), delete(slx); end
new_system(mdl); open_system(mdl);
load_system('sm_lib'); load_system('nesl_utility');

% --- Environment blocks (left column, like the reference report)
hSolv = addb('nesl_utility/Solver Configuration', [mdl '/Solver Configuration'], [40 40 100 80]);
hW    = addb('sm_lib/Frames and Transforms/World Frame', [mdl '/World Frame'], [40 140 100 180]);
hMech = addb('sm_lib/Utilities/Mechanism Configuration', [mdl '/Mechanism Configuration'], [40 240 100 280], ...
    'UniformGravity', 'Constant', 'GravityVector', '[0 -9.80665 0]');

% --- Ground (frame 3เกลอ)
hBase = addSub(mdl, 'Base (Frame)', [190 330 300 390], 'W', '');
buildBase(hBase);

% --- Ground reference frames
hRTc = addRT(mdl, 'Crank Axis O', [190 125 240 165], eye(3), O);
hRTs = addRT(mdl, 'Slider Line',  [190 215 240 255], Ry(90), O + [0; P.e; 0]);

% --- Crank drive (constant RPM): Ramp -> Simulink-PS -> Revolute (input motion)
hRamp = addb('simulink/Sources/Ramp', [mdl '/Crank Angle (deg)'], [190 445 230 475], ...
    'slope', num2str(6*RPM), 'start', '0');
setFirst(hRamp, {'InitialOutput', 'X0'}, num2str(theta0_deg));
hConv = addb('nesl_utility/Simulink-PS Converter', [mdl '/Simulink-PS Converter'], [270 445 300 475], ...
    'Unit', 'deg', 'FilteringAndDerivatives', 'filter', 'SimscapeFilterOrder', '2', ...
    'InputFilterTimeConstant', '1e-3');

hJc = addb('sm_lib/Joints/Revolute Joint', [mdl '/Crank_Joint'], [300 120 360 175], ...
    'TorqueActuationMode', 'ComputedTorque', 'MotionActuationMode', 'InputMotion');

hCrank = addSub(mdl, 'Crank (White)', [410 120 510 175], 'Axis', 'PinB');
buildCrank(hCrank);

hJb = addb('sm_lib/Joints/Revolute Joint', [mdl '/PinB_Joint'], [560 120 620 175], ...
    'PositionTargetSpecify', 'on', 'PositionTargetPriority', 'Low', ...
    'PositionTargetValue', num2str(relB0), 'PositionTargetValueUnits', 'deg');

hRod = addSub(mdl, 'Connecting Rod (Red)', [670 120 770 175], 'PinB', 'PinA');
buildRod(hRod);

hJa = addb('sm_lib/Joints/Revolute Joint', [mdl '/PinA_Joint'], [820 120 880 175]);

hJs = addb('sm_lib/Joints/Prismatic Joint', [mdl '/Slider_Joint'], [300 210 360 265], ...
    'PositionTargetSpecify', 'on', 'PositionTargetPriority', 'Low', ...
    'PositionTargetValue', num2str(x0), 'PositionTargetValueUnits', 'mm');

hSl = addSub(mdl, 'Slider (Pink)', [820 210 920 265], 'Rail', 'PinA');
buildSlider(hSl);

%% ---------------- Connections ----------------
W = pAny(hW, 1);
conn(mdl, W, pAny(hSolv, 1));
conn(mdl, W, pAny(hMech, 1));
conn(mdl, W, pL(hBase, 1));
conn(mdl, W, pL(hRTc, 1));
conn(mdl, W, pL(hRTs, 1));

conn(mdl, pR(hRTc, 1),  pL(hJc, 1));        % O  -> crank joint B
conn(mdl, pR(hJc, 1),   pL(hCrank, 1));     % crank joint F -> crank
conn(mdl, pR(hCrank,1), pL(hJb, 1));        % crank pin -> pin B joint
conn(mdl, pR(hJb, 1),   pL(hRod, 1));
conn(mdl, pR(hRod, 1),  pL(hJa, 1));        % rod end -> pin A joint
conn(mdl, pR(hJa, 1),   pR(hSl, 1));        % pin A joint F -> slider pin
conn(mdl, pR(hRTs, 1),  pL(hJs, 1));        % slider line -> prismatic B
conn(mdl, pR(hJs, 1),   pL(hSl, 1));        % prismatic F -> slider rail frame

phR = get_param(hRamp, 'PortHandles'); phC = get_param(hConv, 'PortHandles');
add_line(mdl, phR.Outport(1), phC.Inport(1), 'autorouting', 'on');
conn(mdl, pR(hConv, 1), pL(hJc, 2));        % motion input q

%% ---------------- Model settings ----------------
Tend = nRev*60/RPM;
safeSet(mdl, 'StopTime', num2str(Tend));
safeSet(mdl, 'Solver', 'ode23t');
safeSet(mdl, 'MaxStep', num2str(min(1e-2, Tend/2000)));
safeSet(mdl, 'RelTol', '1e-6');
safeSet(mdl, 'SimscapeLogType', 'all');
safeSet(mdl, 'SimscapeLogName', 'simlog');
safeSet(mdl, 'ReturnWorkspaceOutputs', 'on');
safeSet(mdl, 'SimMechanicsOpenEditorOnUpdate', 'on');   % เปิด Mechanics Explorer อัตโนมัติ

% Check diagram (auto-fix if the joint's B / q ports come in the other order)
try
    set_param(mdl, 'SimulationCommand', 'update');
catch ME1
    fprintf(2, 'Update failed (%s)\n  -> trying swapped crank-joint ports...\n', ME1.message);
    swapCrankPorts(mdl, hJc, hRTc, hConv);
    try
        set_param(mdl, 'SimulationCommand', 'update');
    catch
        swapCrankPorts(mdl, hJc, hRTc, hConv);  % restore
        rethrow(ME1);
    end
end
save_system(mdl, slx);
fprintf('Saved model: %s\n', slx);

%% ---------------- Run + results ----------------
out = sim(mdl);
plotResults(out, P, RPM, theta0_deg, here);

%% =================== nested builders ===================
    function buildBase(sub)
        pw = pAny(getBlk(sub, 'W'), 1);
        rt = addRT(sub, 'CAD to Model', [150 90 200 130], T.base.R, T.base.t);
        conn(sub, pw, pL(rt, 1));
        parts = {'Frame1_3kelor_smooth', C.frame; 'Frame2_front_face', C.frame2; ...
                 'Frame2_fastener_left_placeholder', C.pin; 'Frame2_fastener_right_placeholder', C.pin};
        if useSTL
            addSolids(sub, rt, parts);
        else  % rail only
            r2 = addRT(sub, 'Rail Center', [150 200 200 240], eye(3), [60; 15; 0]);
            conn(sub, pw, pL(r2, 1));
            b = addb('sm_lib/Body Elements/Brick Solid', [sub '/Rail'], [260 200 310 240], ...
                'BrickDimensions', '[210 10 10]', 'BrickDimensionsUnits', 'mm');
            setSolidColor(b, C.frame);
            conn(sub, pR(r2, 1), pAny(b, 1));
        end
    end

    function buildCrank(sub)
        pa = pAny(getBlk(sub, 'Axis'), 1);  pb = pAny(getBlk(sub, 'PinB'), 1);
        rpin = addRT(sub, 'Crank Radius r', [250 30 300 70], eye(3), [P.r; 0; 0]);
        conn(sub, pa, pL(rpin, 1));  conn(sub, pR(rpin, 1), pb);
        rt = addRT(sub, 'CAD Pose', [150 120 200 160], T.crank.R, T.crank.t);
        conn(sub, pa, pL(rt, 1));
        parts = {'White_crank', C.white; 'Spline_shaft', C.shaft; 'Boss_collar', C.pin; ...
                 'Pin_B_female', C.pin; 'Pin_male_crank', C.pin};
        if useSTL, addSolids(sub, rt, parts); else, addBar(sub, pa, P.r, C.white); end
    end

    function buildRod(sub)
        pa = pAny(getBlk(sub, 'PinB'), 1);  pb = pAny(getBlk(sub, 'PinA'), 1);
        rl = addRT(sub, 'Rod Length L', [250 30 300 70], eye(3), [P.L; 0; 0]);
        conn(sub, pa, pL(rl, 1));  conn(sub, pR(rl, 1), pb);
        rt = addRT(sub, 'CAD Pose', [150 120 200 160], T.rod.R, T.rod.t);
        conn(sub, pa, pL(rt, 1));
        if useSTL, addSolids(sub, rt, {'Red_connecting_link', C.red}); else, addBar(sub, pa, P.L, C.red); end
    end

    function buildSlider(sub)
        pr = pAny(getBlk(sub, 'Rail'), 1);  pp = pAny(getBlk(sub, 'PinA'), 1);
        r2p = addRT(sub, 'Rail to Pin', [150 30 200 70], Ry(-90), [0; 0; 0]);
        conn(sub, pr, pL(r2p, 1));  conn(sub, pR(r2p, 1), pp);
        rt = addRT(sub, 'CAD Pose', [260 120 310 160], T.slider.R, T.slider.t);
        conn(sub, pR(r2p, 1), pL(rt, 1));
        parts = {'Pink_slider', C.pink; 'Pin_A_female', C.pin; 'Pin_male_slider', C.pin};
        if useSTL
            addSolids(sub, rt, parts);
        else
            b = addb('sm_lib/Body Elements/Brick Solid', [sub '/Slider Block'], [380 120 430 160], ...
                'BrickDimensions', '[35 35 35]', 'BrickDimensionsUnits', 'mm');
            setSolidColor(b, C.pink);
            conn(sub, pR(rt, 1), pAny(b, 1));
        end
    end

    function addSolids(sub, rt, parts)
        for k = 1:size(parts, 1)
            y = 40 + 80*(k-1);
            h = addb('sm_lib/Body Elements/File Solid', [sub '/' parts{k,1}], [380 y 440 y+50]);
            safeSet(h, 'ExtGeomFileName', fullfile(partsDir, [parts{k,1} '.stl']));
            safeSet(h, 'UnitType', 'Custom');
            safeSet(h, 'ExtGeomFileUnits', 'mm');
            safeSet(h, 'InertiaType', 'CalculateFromGeometry');
            safeSet(h, 'BasedOnType', 'Density');
            safeSet(h, 'Density', num2str(P.rho));
            safeSet(h, 'DensityUnits', 'kg/m^3');
            setSolidColor(h, parts{k,2});
            conn(sub, pR(rt, 1), pAny(h, 1));
        end
    end

    function addBar(sub, pa, len, col)
        r = addRT(sub, 'Bar Center', [150 220 200 260], eye(3), [len/2; 0; 0]);
        conn(sub, pa, pL(r, 1));
        b = addb('sm_lib/Body Elements/Brick Solid', [sub '/Bar'], [260 220 310 260], ...
            'BrickDimensions', mat2str([len+20 20 8]), 'BrickDimensionsUnits', 'mm');
        setSolidColor(b, col);
        conn(sub, pR(r, 1), pAny(b, 1));
    end
end

%% =================== local helpers ===================
function d = findPartsDir(here)
% หาโฟลเดอร์ Parts_300pct: ข้างไฟล์ .m -> โฟลเดอร์ปัจจุบัน -> ให้ผู้ใช้เลือกเอง
cands = {fullfile(here, 'Parts_300pct'), fullfile(pwd, 'Parts_300pct')};
d = '';
for k = 1:numel(cands)
    if exist(fullfile(cands{k}, 'White_crank.stl'), 'file'), d = cands{k}; return; end
end
if usejava('desktop')
    sel = uigetdir(here, 'เลือกโฟลเดอร์ Parts_300pct (ที่มีไฟล์ .stl)');
    if ischar(sel) && exist(fullfile(sel, 'White_crank.stl'), 'file'), d = sel; end
end
end

function h = addb(src, dst, pos, varargin)
h = add_block(src, dst, 'Position', pos);
for k = 1:2:numel(varargin), safeSet(h, varargin{k}, varargin{k+1}); end
end

function safeSet(h, p, v)
try
    set_param(h, p, v);
catch ME
    if ischar(h), n = h; else, n = getfullname(h); end
    warning('SC:param', 'Could not set "%s" on %s: %s', p, n, ME.message);
end
end

function setSolidColor(h, rgb)
% ตั้งสีให้ Solid block ทุกเวอร์ชันของ MATLAB: หา parameter ที่มีคำว่า Color / VisualProperties เอง
dp = fieldnames(get_param(h, 'DialogParameters'));
ok = false;
for i = 1:numel(dp)
    n = dp{i};
    try
        if ~isempty(regexpi(n, '^(Graphic)?Type$', 'once'))
            set_param(h, n, 'FromGeometry');
        elseif ~isempty(regexpi(n, 'VisualProp', 'once'))
            set_param(h, n, 'Simple');
        end
    catch
    end
end
for i = 1:numel(dp)
    n = dp{i};
    try
        if ~isempty(regexpi(n, 'DiffuseColor$', 'once'))
            set_param(h, n, mat2str(rgb)); ok = true;
        elseif ~isempty(regexpi(n, 'Opacity$', 'once'))
            set_param(h, n, '1');
        end
    catch
    end
end
if ~ok
    warning('SC:color', 'ตั้งสีไม่สำเร็จ: %s', getfullname(h));
end
end

function setFirst(h, names, v)
for k = 1:numel(names)
    try, set_param(h, names{k}, v); return; catch, end
end
warning('SC:param', 'Could not set initial value on %s', getfullname(h));
end

function h = addRT(sys, name, pos, Rm, t)
h = add_block('sm_lib/Frames and Transforms/Rigid Transform', [sys '/' name], 'Position', pos);
Rm = round(Rm, 12); t = round(t(:)', 9);
if isequal(Rm, eye(3))
    safeSet(h, 'RotationMethod', 'None');
else
    safeSet(h, 'RotationMethod', 'RotationMatrix');
    safeSet(h, 'RotationMatrix', mat2str(Rm, 12));
end
if any(t)
    safeSet(h, 'TranslationMethod', 'Cartesian');
    safeSet(h, 'TranslationCartesianOffset', mat2str(t, 10));
    safeSet(h, 'TranslationCartesianOffsetUnits', 'mm');
else
    safeSet(h, 'TranslationMethod', 'None');
end
end

function sub = addSub(mdl, name, pos, leftName, rightName)
sub = [mdl '/' name];
add_block('built-in/Subsystem', sub, 'Position', pos);
hl = add_block('nesl_utility/Connection Port', [sub '/' leftName], 'Position', [30 40 60 60]);
safeSet(hl, 'Side', 'Left');
if ~isempty(rightName)
    hr = add_block('nesl_utility/Connection Port', [sub '/' rightName], 'Position', [560 40 590 60]);
    safeSet(hr, 'Side', 'Right');
end
end

function h = getBlk(sub, name), h = get_param([sub '/' name], 'Handle'); end
function p = pL(h, i), pc = get_param(h, 'PortHandles'); p = pc.LConn(i); end
function p = pR(h, i), pc = get_param(h, 'PortHandles'); p = pc.RConn(i); end
function p = pAny(h, i), pc = get_param(h, 'PortHandles'); a = [pc.LConn pc.RConn]; p = a(i); end
function conn(sys, a, b), add_line(sys, a, b, 'autorouting', 'on'); end

function swapCrankPorts(mdl, hJc, hRTc, hConv)
pc = get_param(hJc, 'PortHandles');
for k = 1:2
    l = get_param(pc.LConn(k), 'Line');
    if l > 0, delete_line(l); end
end
% reconnect the other way round
lb = get_param(pR(hRTc, 1), 'Line'); if lb > 0, delete_line(lb); end
lq = get_param(pR(hConv, 1), 'Line'); if lq > 0, delete_line(lq); end
persistent swapped
if isempty(swapped), swapped = false; end
swapped = ~swapped;
if swapped, iB = 2; iQ = 1; else, iB = 1; iQ = 2; end
conn(mdl, pR(hRTc, 1), pc.LConn(iB));
conn(mdl, pR(hConv, 1), pc.LConn(iQ));
end

function a = wrapTo180local(a), a = mod(a + 180, 360) - 180; end

function plotResults(out, P, RPM, th0, here)
try
    simlog = out.get('simlog');
catch
    simlog = evalin('base', 'simlog');
end
q  = simlog.Crank_Joint.Rz.q.series;
t  = q.time;
th = q.values('deg');
xs = simlog.Slider_Joint.Pz.p.series.values('mm');
vs = simlog.Slider_Joint.Pz.v.series.values('mm/s');
as = gradient(vs, t);

% Analytical slider-crank (offset e), constant omega
w  = RPM*2*pi/60;                        % rad/s
r = P.r; L = P.L; e = P.e;
s  = r*sind(th) - e;  c = cosd(th);
R  = sqrt(L^2 - s.^2);
xa = r*c + R;
dx = -r*sind(th) - (s.*r.*c)./R;                                  % dx/dtheta
d2 = -r*c - ((r*c).^2 - s.*r.*sind(th))./R - (s.*r.*c).^2./R.^3;  % d2x/dtheta2
va = dx*w;   aa = d2*w^2;

f = figure('Name', 'Crank-Slider: Simscape vs Analytical', 'Color', 'w', 'Position', [100 100 900 750]);
subplot(3,1,1); plot(t, xs, 'LineWidth', 1.8); hold on; plot(t, xa, '--', 'LineWidth', 1.4);
ylabel('x_{slider} [mm]'); grid on; legend('Simscape', 'Analytical', 'Location', 'best');
title(sprintf('Crank-Slider (r = %g mm, L = %g mm, e = %g mm), %g RPM, \\theta_0 = %g°', r, L, e, RPM, th0));
subplot(3,1,2); plot(t, vs, 'LineWidth', 1.8); hold on; plot(t, va, '--', 'LineWidth', 1.4);
ylabel('v_{slider} [mm/s]'); grid on;
subplot(3,1,3); plot(t, as, 'LineWidth', 1.8); hold on; plot(t, aa, '--', 'LineWidth', 1.4);
ylabel('a_{slider} [mm/s^2]'); xlabel('time [s]'); grid on;

fprintf('\n===== Results =====\n');
fprintf('Stroke (Simscape)   = %.3f mm  (2r = %.1f mm)\n', max(xs) - min(xs), 2*r);
fprintf('x range             = %.2f ... %.2f mm (from crank axis)\n', min(xs), max(xs));
fprintf('|v| max             = %.1f mm/s\n', max(abs(vs)));
fprintf('|a| max             = %.1f mm/s^2\n', max(abs(as)));
fprintf('max error x         = %.2e mm\n', max(abs(xs - xa)));

try
    Tt = table(t, th, xs, vs, as, xa, va, aa, 'VariableNames', ...
        {'t_s', 'theta_deg', 'x_sim_mm', 'v_sim_mm_s', 'a_sim_mm_s2', 'x_eq_mm', 'v_eq_mm_s', 'a_eq_mm_s2'});
    writetable(Tt, fullfile(here, 'SliderCrank_results.csv'));
    exportgraphics(f, fullfile(here, 'SliderCrank_PVA.png'), 'Resolution', 200);
catch
end
end
