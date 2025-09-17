clear
close all

addpath '0 - Model'\
load("parameters3.mat")


%% Reachability
% Original system
R= ctrb(A,B1);
rank_R= rank(R); % fully reachable if rank 4

% Extended system
Ae= [A zeros(4,1);
     Ce zeros(1,1)];
B1e= [B1 ; De1];
B2e= [B2; De2];

Re= ctrb(Ae,B1e);
rank_Re= rank(Re); % fully reachable if rank 5

%% Observability
% Observability
O= obsv(A,C);
rank_O= rank(O); % fully observable if rank 4
