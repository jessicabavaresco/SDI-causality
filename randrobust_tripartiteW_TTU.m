%   authors:     Jessica Bavaresco
%   last update: May, 2019
%   requires:    Yalmip (https://yalmip.github.io), MOSEK (https://www.mosek.com), and QETLAB (http://www.qetlab.com)

function [q,w_cz_ABC,w_cz_BAC,eta] = randrobust_tripartiteW_TTU(W,M_cz,d)

%   INPUT:    W = tripartite process matrix with Charlie in the future of Alice and Bob
%          M_cz = set of POVMs for Charlie: M_cz(:,:,z,c), where 'z' in the input label and 'c' is the output label
%             d = [dAi dAo dBi dBo dCi], dimension vector, where dAi and dAo are the dimension 
%                 of Alice's input and output spaces, respectively, and analogously for Bob and Charlie.  
%
%   OUTPUT: w_cz_ABC = tripartite process matrix causally ordered from Alice to Bob to Charlie
%           w_cz_BAC = tripartite process matrix causally ordered from Bob to Alice to Charlie
%                  q = parameter for the convex combination of w_cz_ABC and w_cz_BAC, defining a causally separable 
%                      TTU-assemblage w_cz_sep = q w_cz_ABC + (1-q) w_cz_BAC
%                eta = minimum value of the parameter eta such that (1-eta) w_cz + eta I = w_cz_sep

dAi = d(1);
dAo = d(2);
dBi = d(3);
dBo = d(4);
dCi = d(5);
IC  = size(M_cz,3);
OC  = size(M_cz,4);
Id = eye(dAi*dAo*dBi*dBo*dCi)/(dAi*dBi*dCi);

yalmip('clear');

w_cz_ABC = sdpvar(dAi*dAo*dBi*dBo,dAi*dAo*dBi*dBo,IC,OC,'hermitian','complex');
w_cz_BAC = sdpvar(dAi*dAo*dBi*dBo,dAi*dAo*dBi*dBo,IC,OC,'hermitian','complex');
sdpvar eta

F = [];
for z=1:IC
    for c=1:OC
        F = F + [w_cz_ABC(:,:,z,c)>=0,w_cz_BAC(:,:,z,c)>=0];
        F = F + [PartialTrace(kron(eye(dAi*dAo*dBi*dBo),M_cz(:,:,z,c))*((1-eta)*W+eta*Id),2,[dAi*dAo*dBi*dBo dCi])==w_cz_ABC(:,:,z,c)+w_cz_BAC(:,:,z,c)];
    end
    for zz=1:IC
        if zz>z
            F = F + [sum(w_cz_ABC(:,:,z,:),4)==sum(w_cz_ABC(:,:,zz,:),4),sum(w_cz_BAC(:,:,z,:),4)==sum(w_cz_BAC(:,:,zz,:),4)]; 
            F = F + [trace(sum(w_cz_ABC(:,:,z,:),4))==trace(sum(w_cz_ABC(:,:,zz,:),4))];
            F = F + [trace(sum(w_cz_BAC(:,:,z,:),4))==trace(sum(w_cz_BAC(:,:,zz,:),4))];
        end
    end
    F = F + [sum(w_cz_ABC(:,:,z,:),4)==kron(PartialTrace(sum(w_cz_ABC(:,:,z,:),4),2,[dAi*dAo*dBi dBo]),eye(dBo)/dBo),
             PartialTrace(sum(w_cz_ABC(:,:,z,:),4),[3 4],[dAi dAo dBi dBo])==kron(PartialTrace(sum(w_cz_ABC(:,:,z,:),4),[2 3 4],[dAi dAo dBi dBo]),eye(dAo)/dAo)];
    F = F + [sum(w_cz_BAC(:,:,z,:),4)==PermuteSystems(kron(PartialTrace(sum(w_cz_BAC(:,:,z,:),4),2,[dAi dAo dBi dBo]),eye(dAo)/dAo),[1 4 2 3],[dAi dBi dBo dAo]),
             PartialTrace(sum(w_cz_BAC(:,:,z,:),4),[1 2],[dAi dAo dBi dBo])==kron(PartialTrace(sum(w_cz_BAC(:,:,z,:),4),[1 2 4],[dAi dAo dBi dBo]),eye(dBo)/dBo)];
end

flag = solvesdp(F,eta,sdpsettings('solver','mosek','verbose',0))


w_cz_ABC = double(w_cz_ABC);
w_cz_BAC = double(w_cz_BAC);
eta      = double(eta);

q = trace(sum(w_cz_ABC(:,:,1,:),4))/(dAo*dBo);
for z=1:IC
    for c=1:OC
        w_cz_ABC(:,:,z,c) = w_cz_ABC(:,:,z,c)/q;
        w_cz_BAC(:,:,z,c) = w_cz_BAC(:,:,z,c)/(1-q);
    end
end

        


