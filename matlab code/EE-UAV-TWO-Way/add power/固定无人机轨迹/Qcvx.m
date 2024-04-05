function [q,distance2] = Qcvx(theta,E_VU,E_VR,W,Y,L_VR,rho_l,omega_l,distance_l_SCA,CARposition,...
                          VUAV_max,det,T,M,roadlength,roadwidth,UAVposition) 
% cvx_clear
% %cvx_solver SeDuMi
% cvx_begin
%             variable q(T,2)
q=UAVposition;
            for t=1:T
                for m=1:M
                    distance2(t,m)=square_pos(norm([q(t,1) q(t,2)]-[CARposition(m,1,t) CARposition(m,2,t)]));
                end
            end
            L_A=sum(sum(W*Y'.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%凸+凹？
            L_total=L_A+L_VR;
 %           minimize (-1*L_A)
%              minimize (-0.0001*(L_total/(theta*E_VR+1*E_VU)))   %(1-theta)
%              /(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl));%原来的
                
%             subject to
%      
%             for t=1:T-1
%                 norm([q(t,1) q(t,2)]-[q(t+1,1) q(t+1,2)]) <= VUAV_max*det;
%             end
            q(1,1)==0;
            q(1,2)==0;
            q(T,1)==roadlength;
            q(T,2)==roadwidth;%起点终点的约束

% cvx_end                       
                      