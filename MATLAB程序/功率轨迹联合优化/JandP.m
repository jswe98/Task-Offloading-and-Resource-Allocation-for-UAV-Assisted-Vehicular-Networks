syms P1 P2 P3 Pmax
syms a b c q p
syms x
syms L
syms c2
N0=96
gama=1e-12
B=60
Pmax=26
P1(1)=26/3
P2(1)=26/3
P3(1)=Pmax-P1(1)-P2(1)
xx(1)=400
dmu(1)=((B-xx(1)).^2+100.^2).^(1./2)
dud(1)=(xx(1).^2+100.^2).^(1./2)
Pout(1)=eval(1-exp(-N0.*gama./2.*(20./(P3(1).*5.^(-4))+1./(P1(1).*dmu(1).^(-4))+1./(P2(1).*dud(1).^(-4)))))
for t=1:5
   syms lambda
   a=4./P1(t)+4./P2(t)
   b=-12.*B./P1(t)
   c=4.*(3.*B.^2+100.^2)./P1(t)+4.*100.^2./P2(t)
   q=-lambda./a+2.*b.^3./(27.*a.^3)-b.*c./(3.*a.^2)
   p=c./a-b.^2./(3.*a.^2)
   x(t)=(-q./2+((q./2).^2+(p./3).^3).^(1./2)).^(1./3)+(-q./2-((q./2).^2+(p./3).^3).^(1./2)).^(1./3)-b./(3.*a)
   L(t)=(5.^4./((Pmax-P1(t)-P2(t))./20)).*20+((B-x(t)).^2+100.^2).^2./P1(t)+(x(t).^2+100.^2)./P2(t)+lambda.*(xx(t)-x(t)-0.05)
   y(t)=diff(L(t))
   lambda=50000
   alpha=100
   kesi=200
   num=0
   while 1
    num=num+1
    tidu(t)=abs(eval(subs(y(t),lambda)))
    if tidu(t)>kesi
        lambda=lambda+alpha*tidu(t)
    end
    if tidu(t)<=kesi
        break
    end
   end
   xx(t+1)=abs(eval(subs(x(t),lambda)))
   c2=((B-xx(t+1)).^2+100.^2)./(xx(t+1).^2+100.^2)
   P1(t+1)=vpa((c2./(1+c2)).*Pmax)
   P2(t+1)=vpa((0.8./(1+c2)).*Pmax)
   P3(t+1)=vpa((Pmax-P1(t+1)-P2(t+1))./20)
   dmu(t+1)=((B-xx(t+1)).^2+100.^2).^(1./2)
   dud(t+1)=(xx(t+1).^2+100.^2).^(1./2)
   Pout(t+1)=eval(1-exp(-N0.*gama./2.*(20./(P3(t+1).*5.^(-4))+1./(P1(t+1).*dmu(t+1).^(-4))+1./(P2(t+1).*dud(t+1).^(-4)))))
   t=t+1
   clear lambda
end
t=0:5
i=0:100:500
figure(1)
plot(i,xx(t+1),'blue-v','linewidth',1.5)
xlabel('\fontname{宋体}时隙')
ylabel('\fontname{宋体}无人机到基站的距离\fontname{Times New Roman}(m)')
legend('\fontname{宋体}无人机轨迹')
set(gca,'FontSize',13)
grid on
figure(2)
plot(i,P1(t+1),'red-v','linewidth',1.5)
hold on
plot(i,P2(t+1),'blue-p','linewidth',1.5)
hold on
plot(i,P3(t+1),'green-o','linewidth',1.5)
xlabel('\fontname{宋体}时隙')
ylabel('\fontname{宋体}发射功率\fontname{Times New Roman}(dBm)')
legend('\fontname{宋体}簇头车辆发射功率','\fontname{宋体}无人机发射功率','\fontname{宋体}簇内车辆发射功率')
set(gca,'FontName','Times New Roman','FontSize',13)
grid on
axis([0 500 -1 20])