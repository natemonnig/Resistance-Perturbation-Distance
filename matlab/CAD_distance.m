function [ d ] = CAD_distance(A0,A1,R0,R1)
% Commute Time Anomaly Detection distance from
% Kumar Sricharan and Kamalika Das. 2014. Localizing anomalous changes in time­evolving graphs. In 1283 Proceedings of the 2014 ACM International Conference on Management of Data (SIGMOD). ACM, 1347­ 1284 1358.
% The distance computed here is for an empty edge set S (which measures the
% graph change due to changes in the full edge set)

commute_0 = sum(sum(A0))*R0;
commute_1 = sum(sum(A1))*R1;
d=sum(sum( abs(A0-A1).*abs(commute_0 - commute_1) ));

end

