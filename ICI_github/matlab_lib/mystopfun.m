function stopnow = mystopfun(problem, x, info, last)
    stopnow = (last >= 3 && info(last-2).cost - info(last).cost < 1e-6);
%     stopnow = (info(last-2).cost - info(last).cost < 1e-3);
end