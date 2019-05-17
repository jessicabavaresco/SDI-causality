%   authors:     Jessica Bavaresco
%   last update: May, 2019
%   requires:    Yalmip (https://yalmip.github.io) and QETLAB (http://www.qetlab.com)

function [q,W_ABC,W_BAC,eta] = randrobust_tripartiteW_TTT(W,d)

%   INPUT: W = tripartite process matrix with Charlie in the future of Alice and Bob
%          d = [dAi dAo dBi dBo dCi], dimension vector, where dAi and dAo are the dimension 
%              of Alice's input and output spaces, respectively, and analogously for Bob and Charlie.  
%
%   OUTPUT: W_ABC = tripartite process matrix causally ordered from Alice to Bob to Charlie
%           W_BAC = tripartite process matrix causally ordered from Bob to Alice to Charlie
%               q = parameter for the convex combination of W_ABC and W_ABC, defining a causally 
%                   separable tripartite process matrix Wsep = q W_ABC + (1-q) W_BAC
%             eta = minimum value of the parameter eta such that (1-eta) W + eta I = Wsep

dAi = d(1);
dAo = d(2);
dBi = d(3);
dBo = d(4);
dCi = d(5);
Id = eye(dAi*dAo*dBi*dBo*dCi)/(dAi*dBi*dCi);

yalmip('clear');

W_ABC = sdpvar(dAi*dAo*dBi*dBo*dCi,dAi*dAo*dBi*dBo*dCi,'hermitian','complex');
W_BAC = sdpvar(dAi*dAo*dBi*dBo*dCi,dAi*dAo*dBi*dBo*dCi,'hermitian','complex');
sdpvar eta

F = [W_ABC>=0,W_BAC>=0,
    PartialTrace(W_ABC,5,[dAi dAo dBi dBo dCi])==kron(PartialTrace(W_ABC,[4 5],[dAi dAo dBi dBo dCi]),eye(dBo)/dBo),
    PartialTrace(W_ABC,[3 4 5],[dAi dAo dBi dBo dCi])==kron(PartialTrace(W_ABC,[2 3 4 5],[dAi dAo dBi dBo dCi]),eye(dAo)/dAo),
    PartialTrace(W_BAC,5,[dAi dAo dBi dBo dCi])==PermuteSystems(kron(PartialTrace(W_BAC,[2 5],[dAi dAo dBi dBo dCi]),eye(dAo)/dAo),[1 4 2 3],[dAi dBi dBo dAo]),
    PartialTrace(W_BAC,[1 2 5],[dAi dAo dBi dBo dCi])==kron(PartialTrace(W_BAC,[1 2 4 5],[dAi dAo dBi dBo dCi]),eye(dBo)/dBo),
    (1-eta)*W+eta*Id==W_ABC+W_BAC];
    
flag_robustSWITCH = solvesdp(F,eta,sdpsettings('solver','mosek','verbose',0))

eta = double(eta);
W_ABC  = double(W_ABC);
W_BAC  = double(W_BAC);

q = trace(W_ABC);
W_ABC = dAo*dBo*W_ABC/trace(W_ABC);
W_BAC = dAo*dBo*W_BAC/trace(W_BAC);


