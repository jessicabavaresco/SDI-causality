%   authors:     Jessica Bavaresco
%   last update: May, 2019
%   requires:    Yalmip (https://yalmip.github.io) and QETLAB (http://www.qetlab.com)

function [eta_TTT,eta_UTT,eta_TTU,eta_TUU] = switch_randrobustness

ket_0     = [1;0];
ket_1     = [0;1];
ket_plus  = (1/sqrt(2))*(ket_1+ket_0);
ket_minus = (1/sqrt(2))*(-ket_1+ket_0);

d = [2 2 2 2 2];

%switch process' parameters
psi   = ket_0;
alpha = 1/sqrt(2);
beta  = 1/sqrt(2);

phi    = sqrt(2)*MaxEntangled(2);
SWITCH = alpha*Tensor(psi,phi,phi,ket_0) + beta*PermuteSystems(Tensor(psi,phi,phi,ket_1),[2 1 3],[4 4 4]);
SWITCH = SWITCH*SWITCH';

%reduced switch process for given psi, alpha, and beta
W = PartialTrace(SWITCH,2,[16 2 2]);

%Alice's instruments
A_ax(:,:,1,1) = kron(ket_0*ket_0',ket_0*ket_0');
A_ax(:,:,1,2) = kron(ket_1*ket_1',ket_1*ket_1');
A_ax(:,:,2,1) = kron(ket_plus*ket_plus',ket_plus*ket_plus');
A_ax(:,:,2,2) = kron(ket_minus*ket_minus',ket_minus*ket_minus');

%Bob's instruments
B_by = A_ax;

%Charlie's POVM
M_cz(:,:,1,1) = ket_plus*ket_plus';
M_cz(:,:,1,2) = ket_minus*ket_minus';


% calculates the bounds for \eta^* presented on Table II:

%TTT
[q_TTT,W_ABC,W_BAC,eta_TTT]           = randrobust_tripartiteW_TTT(W,d);                      %#ok<*ASGLU>
%UTT
[q_UTT,w_ax_ABC,w_ax_BAC,eta_UTT]     = randrobust_tripartiteW_UTT(W,A_ax,d);
%TTU
[q_TTU,w_cz_ABC,w_cz_BAC,eta_TTU]     = randrobust_tripartiteW_TTU(W,M_cz,d);
%TUU
[q_TUU,w_bcyz_ABC,w_bcyz_BAC,eta_TUU] = randrobust_tripartiteW_TUU(W,B_by,M_cz,d);
