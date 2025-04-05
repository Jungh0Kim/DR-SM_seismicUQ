% Jungho Kim
% Set of ground motion features used as physics-based input variables (Table 1)

function [PDR_vars_sv] = PDR_GM(GM_cell)

g = 9.81;     % m/s^2
n_sim = size(GM_cell,1);

PDR_vars_sv = [];
for i_sim = 1:n_sim
    GM_sv = GM_cell{i_sim,1};
    t_hist = GM_sv.t_hist;

    GM = GM_sv.GM;  
    dt = t_hist(10) - t_hist(9);
    GV = cumtrapz(t_hist, GM.*g);
    GD = cumtrapz(t_hist, GV);

    PGA = max(abs(GM));
    PGV = max(abs(GV));
    PGD = max(abs(GD));

    sPeriod = [0.01,0.02,0.022,0.025,0.029,0.03,0.032,0.035,0.036,...
        0.04,0.042,0.044,0.045,0.046,0.048,0.05,0.055,0.06,0.065,0.067,0.07,...
        0.075,0.08,0.085,0.09,0.095,0.1,0.11,0.12,0.125,0.13,0.133,0.14,0.15,...
        0.16,0.17,0.18,0.19,0.2,0.22,0.24,0.25,0.26,0.28,0.29,0.3,0.32,0.34,...
        0.35,0.36,0.38,0.4,0.42,0.44,0.45,0.46,0.48,0.5,0.55,0.6,0.65,0.667,...
        0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,...
        2,2.2,2.4,2.5,2.6,2.8,3,3.2,3.4,3.5,3.6,3.8,4,4.2,4.4,4.6,4.8,5,7.5,10];
    xi = 0.05;  % critical damping
    [Pseudo_SA, Pseudo_SV, SD, SA, SV] = responseSpectra(xi, sPeriod, GM, dt);
    SA = SA./g;            % g
    
    AriasTT_hist = cumtrapz(t_hist, (GM.*g).^2).*pi./(2*g);
    Arias_TT = AriasTT_hist(end);
    AriasTT_N_hist = AriasTT_hist./Arias_TT;
    Arias_TT_mean = mean(AriasTT_hist);
    Arias_TT_median = median(AriasTT_hist);

    % DURATION
    time_dur1 = t_hist(AriasTT_hist >= 0.05*Arias_TT & AriasTT_hist <= 0.75*Arias_TT);
    t_5_75 = [time_dur1(1), time_dur1(end)];
    D_5_75 = time_dur1(end) - time_dur1(1);
    time_dur2 = t_hist(AriasTT_hist >= 0.05*Arias_TT & AriasTT_hist <= 0.95*Arias_TT);
    t_5_95 = [time_dur2(1), time_dur2(end)];
    D_5_95 = time_dur2(end) - time_dur2(1);
    time_dur45 = t_hist(AriasTT_hist >= 0.05*Arias_TT & AriasTT_hist <= 0.45*Arias_TT);
    t_45 = [time_dur45(end)];

    [~, I_mean] = min(abs(AriasTT_hist - Arias_TT_mean));
    t_Arias_mean = t_hist(I_mean,1);
    [~, I_median] = min(abs(AriasTT_hist - Arias_TT_median));
    t_Arias_median = t_hist(I_median,1);

    PDR_vars = [PGA, PGV, PGD, SA, Arias_TT, Arias_TT_mean, Arias_TT_median, D_5_75, D_5_95, t_45, t_Arias_mean, t_Arias_median];
    PDR_vars_sv = [PDR_vars_sv; PDR_vars];
end

end % function end
