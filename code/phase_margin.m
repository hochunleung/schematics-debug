pkg load control;

% 電路參數
gm1 = 2.5e-3; R1 = 8e3; C1 = 100e-15; 
gm2 = 10e-3; R2 = 5e3; C2 = 2e-12;   
Cc = 100e-15; 
Rz = 200;
Cgd2 = 0;
beta = 1;
A0 = (gm1 * R1) * (gm2 * R2);
s = tf('s');

% 自動計算 Rz (抵消次極點)
b_tmp = R1*(C1 + Cc) + R2*(C2 + Cc) + gm2*R2*R1*Cc;
a_tmp = R1*R2*(C1*C2 + Cc*C1 + Cc*C2);
p_roots = roots([a_tmp, b_tmp, 1]);
p2_approx = max(abs(p_roots));
Rz_cal = (1 / (p2_approx * Cc)) + (1 / gm2);

% 三階係數
%b = R1*(C1 + Cc) + R2*(C2 + Cc) + gm2*R2*R1*Cc;
%a = R1*R2*(C1*C2 + Cc*C1 + Cc*C2) + Rz*Cc*(R1*C1 + R2*C2);
%c = Rz*R1*R2*C1*C2*Cc;

Ctot = Cc + Cgd2;
c = R1*(C1 + Ctot) + R2*(C2 + Ctot) + gm2*R2*R1*Ctot;
b = R1*R2*(C1*C2 + Ctot*C1 + Ctot*C2) + Rz*Cc*(R1*C1 + R2*C2); % b項中，只有Cc與Rz耦合
a = Rz*R1*R2*C1*C2*Ctot; % Ctot?

% 求解三個極點
all_poles = roots([a, b, c, 1]);
f_poles = sort(abs(all_poles)) / (2*pi); % 轉換為 Hz 並排序

fprintf('\n--- 三階模型極點詳細分析 --- \n');
fprintf('Cc = %.2f fF, Rz = %.2f Ohm\n', Cc*1e15, Rz);
fprintf('主極點   f_p1: %.2f MHz\n', f_poles(1)/1e6);
fprintf('次極點   f_p2: %.2f MHz\n', f_poles(2)/1e6);
fprintf('第三極點 f_p3: %.2f MHz\n', f_poles(3)/1e6);

% 零點
z1 = 1 / (Cc * (1/gm2 - Rz) + Cgd2/gm2);
fprintf('LHP 零點 f_z : %.2f MHz\n', abs(z1)/(2*pi*1e6));

% 穩定性
H = A0 * (1 - s/z1) / (a*s^3 + b*s^2 + c*s + 1);
T = beta * H;
[Gm, Pm, Wcg, Wgc] = margin(T);
G_close = feedback(H, beta);
fprintf('環路截止頻率 (f_gc): %.2f MHz\n', Wgc / (2*pi*1e6));
fprintf('相位裕度 PM  : %.2f 度\n', Pm);
fprintf('增益裕度 PM  : %.2f dB\n', Gm);

% 繪圖對比
figure;
bode(T); title('環路傳遞函數 T(s) 的 Bode 圖 (用於分析穩定性)');
grid on;

%figure;
%step(G_close); title('閉環系統 G_{CL}(s) 的階躍響應');
%grid on;