EEno_f=[5.45240732378429;5.19454125049931;5.15307851286512;5.14683997913937;5.14605649148186;5.14607964273525;5.14620003517403;5.14631622631429;5.14641787042807;5.14650681610722;5.14658566695020;5.14665645104656;5.14672066464466;5.14677942615514;5.14683359086836;5.14688382717413;5.14693066733951;5.14697454216713;5.14701580527604;5.14705475052161;5.14709162475804;5.14712663736360;5.14715996746709;5.14719176951126;5.14722217759257;5.14725130888689;5.14727926638285;5.14730614108452;5.14733201380229;5.14735695662118;5.14738103411385;5.14740430434980;5.14742681974061;5.14744862775206;5.14746977150763;5.14749029030268;5.14751022004461;5.14752959363161;5.14754844127980;5.14756679080720;5.14758466788097;5.14760209623363;5.14761909785280;5.14763569314817;5.14765190109900;5.14766773938477;5.14768322450122;5.14769837186378;5.14771319589985;5.14772771013155];
EE=[7.92530418326120;8.01307398464635;8.03034039642195;8.03407850798607;8.03510640090783;8.03554130090161;8.03581584977035;8.03602801822516;8.03620551354404;8.03635929761186;8.03649525479467;8.03661716530881;8.03672768402350;8.03682876812015;8.03692190715223;8.03700826244528;8.03708875738922;8.03716413846851;8.03723501783303;8.03730190376047;8.03736522293226;8.03742533703340;8.03748255533022;8.03753714434278;8.03758933538260;8.03763933049770;8.03768730721267;8.03773342234570;8.03777781511032;8.03782060965699;8.03786191717171;8.03790183762134;8.03794046121459;8.03797786963260;8.03801413707141;8.03804933112976;8.03808351356905;8.03811674096684;8.03814906528141;8.03818053434147;8.03821119227258;8.03824107986999;8.03827023492562;8.03829869251588;8.03832648525585;8.03835364352432;8.03838019566375;8.03840616815832;8.03843158579296;8.03845647179568];
EEno_p=[0;7.79034301527730;8.64397651463452;4.30018270244318;0.950428727989293;1.51100596202124;3.13587980677853;4.46090225248120;5.23207727519971;5.54518996398638;5.59627763066812;5.55149152035170;5.49703330628648;5.45822418931528;5.43472259839799;5.42098406242914;5.41270662927856;5.40742995612865;5.40386946785540;5.40135231366673;5.39950784729766;5.39811832311858;5.39704807627101;5.39620849694101;5.39553953821576;5.39499929776680;5.39455782021900;5.39419324794583;5.39388934672037;5.39363386947591;5.39341744800862;5.39323282648896;5.39307432154786;5.39293743568497;5.39281857634936;5.39271484906094;5.39262390318831;5.39254381568706;5.39247300254661;5.39241015069567;5.39235416517399;5.39230412780856;5.39225926463898;5.39221892005383;5.39218253611589;5.39214963592975;5.39211981018066;5.39209270617825;5.39206801889095;5.39204548357194];
figure
plot((EE(:,1)),'-+r');
xlabel('iteration');
ylabel('EE');
hold on
plot((EEno_f(:,1)),'-+b');
xlabel('iteration');
ylabel('EE');
plot((EEno_p(:,1)),'-+g');
legend('EE','EEno_f','EEno_p')
xlabel('iteration');
set(gca,'ylim',[0,12]); 
ylabel('EE');