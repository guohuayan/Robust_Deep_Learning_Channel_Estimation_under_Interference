function [phase] = eliminate_delay_cluster(phase)
    thr=-25;
%     thr=-50;
    thr=10^(thr/10);
    p_phase=abs(phase).^2;
    m_phase=max(p_phase)*thr;
    idx=(p_phase>=m_phase);
    phase=phase.*idx;
end