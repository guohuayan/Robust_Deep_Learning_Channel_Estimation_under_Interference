function [support_w] = eliminate_delay_cluster_group(support_w)
    thr=-25;
%     thr=-50;
    thr=10^(thr/10);
    m_phase=max(support_w)*thr;
    idx=(support_w>=m_phase);
    support_w=support_w.*idx;
end