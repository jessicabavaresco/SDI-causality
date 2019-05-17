%   authors:     Jessica Bavaresco
%   last update: May, 2019
%   requires:    Yalmip (https://yalmip.github.io), MOSEK (https://www.mosek.com), and QETLAB (http://www.qetlab.com)

function [q,w_ax_ABC,w_ax_BAC,eta] = randrobust_tripartiteW_UTT(W,A_ax,d)

%   INPUT:    W = tripartite process matrix with Charlie in the future of Alice and Bob
%          A_ax = set of instruments for Alice: A_ax(:,:,x,a), where 'x' in the input label and 'a' is the output label
%             d = [dAi dAo dBi dBo dCi], dimension vector, where dAi and dAo are the dimension 
%                 of Alice's input and output spaces, respectively, and analogously for Bob and Charlie.  
%
%   OUTPUT: w_ax_ABC = tripartite process matrix causally ordered from Alice to Bob to Charlie
%           w_ax_BAC = tripartite process matrix causally ordered from Bob to Alice to Charlie
%                  q = parameter for the convex combination of w_ax_ABC and w_ax_BAC, defining a causally separable 
%                      TTU-assemblage w_ax_sep = q w_ax_ABC + (1-q) w_ax_BAC
%                eta = minimum value of the parameter eta such that (1-eta) w_ax + eta I = w_ax_sep

dAi = d(1);
dAo = d(2);
dBi = d(3);
dBo = d(4);
dCi = d(5);
IA  = size(A_ax,3);
OA  = size(A_ax,4);
Id = eye(dAi*dAo*dBi*dBo*dCi)/(dAi*dBi*dCi);

yalmip('clear');

w_ax_ABC = sdpvar(dBi*dBo*dCi,dBi*dBo*dCi,IA,OA,'hermitian','complex');
w_ax_BAC = sdpvar(dBi*dBo*dCi,dBi*dBo*dCi,IA,OA,'hermitian','complex');
sdpvar eta

F = [];
for x=1:IA
    F = F + [PartialTrace(sum(w_ax_BAC(:,:,x,:),4),3,[dBi dBo dCi])==kron(PartialTrace(sum(w_ax_BAC(:,:,x,:),4),[2 3],[dBi dBo dCi]),eye(dBo)/dBo)];
    for a=1:OA
        F = F + [w_ax_ABC(:,:,x,a)>=0,w_ax_BAC(:,:,x,a)>=0];
        F = F + [PartialTrace(w_ax_ABC(:,:,x,a),3,[dBi dBo dCi])==kron(PartialTrace(w_ax_ABC(:,:,x,a),[2 3],[dBi dBo dCi]),eye(dBo)/dBo)];       
        F = F + [PartialTrace(Tensor(A_ax(:,:,x,a),eye(dBi*dBo),eye(dCi))*((1-eta)*W+eta*Id),1,[dAi*dAo dBi*dBo dCi])==w_ax_ABC(:,:,x,a)+w_ax_BAC(:,:,x,a)];
    end
    for xx=1:IA
        if xx>x
            F = F + [PartialTrace(sum(w_ax_BAC(:,:,x,:),4),3,[dBi dBo dCi])==PartialTrace(sum(w_ax_BAC(:,:,xx,:),4),3,[dBi dBo dCi])];
            F = F + [trace(sum(w_ax_ABC(:,:,x,:),4))==trace(sum(w_ax_ABC(:,:,xx,:),4))];
            F = F + [trace(sum(w_ax_BAC(:,:,x,:),4))==trace(sum(w_ax_BAC(:,:,xx,:),4))];
        end
    end
end

flag = solvesdp(F,eta,sdpsettings('solver','mosek','verbose',0))

w_ax_ABC = double(w_ax_ABC);
w_ax_BAC = double(w_ax_BAC);
eta      = double(eta);

q = trace(sum(w_ax_ABC(:,:,1,:),4))/dBo;
for x=1:IA
    for a=1:OA
        w_ax_ABC(:,:,x,a) = w_ax_ABC(:,:,x,a)/q;
        w_ax_BAC(:,:,x,a) = w_ax_BAC(:,:,x,a)/(1-q);
    end
end
