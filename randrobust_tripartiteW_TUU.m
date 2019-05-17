%   authors:     Jessica Bavaresco
%   last update: May, 2019
%   requires:    Yalmip (https://yalmip.github.io), MOSEK (https://www.mosek.com), and QETLAB (http://www.qetlab.com)

function [q,w_bcyz_ABC,w_bcyz_BAC,eta] = randrobust_tripartiteW_TUU(W,B_by,M_cz,d)

%   INPUT:    W = tripartite process matrix with Charlie in the future of Alice and Bob
%          B_by = set of instruments for bob: B_by(:,:,y,b), where 'y' in the input label and 'b' is the output label
%          M_cz = set of POVMs for Charlie: M_cz(:,:,z,c), where 'z' in the input label and 'c' is the output label
%             d = [dAi dAo dBi dBo dCi], dimension vector, where dAi and dAo are the dimension 
%                 of Alice's input and output spaces, respectively, and analogously for Bob and Charlie.  
%
%   OUTPUT: w_bcyz_ABC = tripartite process matrix causally ordered from Alice to Bob to Charlie
%           w_bcyz_BAC = tripartite process matrix causally ordered from Bob to Alice to Charlie
%                    q = parameter for the convex combination of w_cz_ABC and w_cz_BAC, defining a causally separable 
%                        TTU-assemblage w_bcyz_sep = q w_bcyz_ABC + (1-q) w_bcyz_BAC
%                   eta = minimum value of the parameter eta such that (1-eta) w_bcyz + eta I = w_bcyz_sep

dAi = d(1);
dAo = d(2);
dBi = d(3);
dBo = d(4);
dCi = d(5);
IB  = size(B_by,3);
OB  = size(B_by,4);
IC  = size(M_cz,3);
OC  = size(M_cz,4);
Id  = eye(dAi*dAo*dBi*dBo*dCi)/(dAi*dBi*dCi);

yalmip('clear');

w_bcyz_ABC = sdpvar(dAi*dAo,dAi*dAo,IB,IC,OB,OC,'hermitian','complex');
w_bcyz_BAC = sdpvar(dAi*dAo,dAi*dAo,IB,IC,OB,OC,'hermitian','complex');
sdpvar eta

F = [];

for y=1:IB
    for z=1:IC
        for b=1:OB
            for c=1:OC
                F = F + [w_bcyz_ABC(:,:,y,z,b,c)>=0,w_bcyz_BAC(:,:,y,z,b,c)>=0];
                F = F + [PartialTrace(Tensor(eye(dAi*dAo),B_by(:,:,y,b),M_cz(:,:,z,c))*((1-eta)*W+eta*Id),[2 3],[dAi*dAo dBi*dBo dCi])==w_bcyz_ABC(:,:,y,z,b,c)+w_bcyz_BAC(:,:,y,z,b,c)];
            end 
        end
        F = F + [sum(sum(w_bcyz_ABC(:,:,y,z,:,:),6),5)==kron(PartialTrace(sum(sum(w_bcyz_ABC(:,:,y,z,:,:),6),5),2,[dAi dAo]),eye(dAo)/dAo)];
    end
end


for y=1:IB
    for b=1:OB
        for z=1:IC
            for zz=1:IC
                if zz>z
                    F = F + [sum(w_bcyz_ABC(:,:,y,z,b,:),6)==sum(w_bcyz_ABC(:,:,y,zz,b,:),6),sum(w_bcyz_BAC(:,:,y,z,b,:),6)==sum(w_bcyz_BAC(:,:,y,zz,b,:),6)];
                end
            end
            F = F + [sum(w_bcyz_BAC(:,:,y,z,b,:),6)==kron(PartialTrace(sum(w_bcyz_BAC(:,:,y,z,b,:),6),2,[dAi dAo]),eye(dAo)/dAo)];
        end
    end
end
     
for y=1:IB
    for z=1:IC
        for yy=1:IB
           for zz=1:IC
               if (yy-1)*IC+zz>(y-1)*IC+z
                  F = F + [sum(sum(w_bcyz_ABC(:,:,y,z,:,:),6),5)==sum(sum(w_bcyz_ABC(:,:,yy,zz,:,:),6),5)];
                  F = F + [trace(sum(sum(w_bcyz_ABC(:,:,y,z,:,:),6),5))==trace(sum(sum(w_bcyz_ABC(:,:,yy,zz,:,:),6),5))];
                  F = F + [trace(sum(sum(w_bcyz_BAC(:,:,y,z,:,:),6),5))==trace(sum(sum(w_bcyz_BAC(:,:,yy,zz,:,:),6),5))];
               end
           end
        end
    end
end

flag = solvesdp(F,eta,sdpsettings('solver','mosek','verbose',0))

w_bcyz_ABC = double(w_bcyz_ABC);
w_bcyz_BAC = double(w_bcyz_BAC);
eta        = double(eta);

q = trace(sum(sum(w_bcyz_ABC(:,:,1,1,:,:),6),5))/dAo;
for y=1:IB
    for z=1:IC
        for b=1:OB
            for c=1:OC
                w_bcyz_ABC(:,:,y,z,b,c) = w_bcyz_ABC(:,:,y,z,b,c)/q;
                w_bcyz_BAC(:,:,y,z,b,c) = w_bcyz_BAC(:,:,y,z,b,c)/(1-q);
            end
        end
    end
end
