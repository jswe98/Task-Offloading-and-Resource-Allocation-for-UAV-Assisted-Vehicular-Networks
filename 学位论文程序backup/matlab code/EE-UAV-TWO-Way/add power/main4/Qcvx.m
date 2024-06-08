function [q,distance2] = Qcvx(theta,E_VU,E_VR,W,Y,L_VR,rho_l,omega_l,distance_l_SCA,CARposition,...
                          VUAV_max,det,T,M) 
cvx_clear
%cvx_solver SeDuMi
cvx_begin
            variable q(T,2)
            for t=1:T
                for m=1:M
                    distance2(t,m)=square_pos(norm([q(t,1) q(t,2)]-[CARposition(m,1,t) CARposition(m,2,t)]));
                end
            end
            L_A=sum(sum(W*Y'.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%Í¹+°¼£¿
            L_total=L_A+L_VR;
 %           minimize (-1*L_A)
             minimize (-0.0001*(L_total/(theta*E_VR+(1-theta)*E_VU)))
%              /(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl));%Ô­À´µÄ
                
            subject to
     
            for t=1:T-1
                norm([q(t,1) q(t,2)]-[q(t+1,1) q(t+1,2)]) <= VUAV_max*det;
            end
            q(1,1)==0;
            q(1,2)==0;

cvx_end                       
                      