figure;
plot(rhoSample,neSample,'b');
hold on;
plot(rhoSample,polyval(fliplr(nepoly), rhoSample),'r');
legend('original ne', 'fitted ne'); 
title('ne')
figure;
plot(rhoSample,polyval(polyder(fliplr(nepoly)), rhoSample),'r');
title('fitted dne/drho')

if IonType~=3
    figure;
    plot(rhoSample,niSample,'b');
    hold on;
    plot(rhoSample,polyval(fliplr(nipoly), rhoSample),'r');
    legend('original ni', 'fitted ni'); 
    title('ni')
    figure;
    plot(rhoSample,polyval(polyder(fliplr(nipoly)), rhoSample),'r');
    title('fitted dni/drho')
    if IonType==1
        figure;
        plot(rhoSample,TiSample,'b');
        hold on;
        plot(rhoSample,polyval(fliplr(Tipoly), rhoSample),'r');
        legend('original Ti', 'fitted Ti'); 
        title('Ti')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Tipoly)), rhoSample),'r');
        title('fitted dTi/drho')
    elseif IonType==2
        figure;
        plot(rhoSample,PiSample,'b');
        hold on;
        plot(rhoSample,polyval(fliplr(Pipoly), rhoSample),'r');
        legend('original Pi', 'fitted Pi'); 
        title('Pi')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Pipoly)), rhoSample),'r');
        title('fitted dPi/drho')
    end
end


if AlphaType~=3
    figure;
    plot(rhoSample,naSample,'b');
    hold on;
    plot(rhoSample,polyval(fliplr(napoly), rhoSample),'r');
    legend('original na', 'fitted na'); 
    title('na')
    figure;
    plot(rhoSample,polyval(polyder(fliplr(napoly)), rhoSample),'r');
    title('fitted dna/drho')
    if AlphaType==1
        figure;
        plot(rhoSample,TaSample,'b');
        hold on;
        plot(rhoSample,polyval(fliplr(Tapoly), rhoSample),'r');
        legend('original Ta', 'fitted Ta'); 
        title('Ta')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Tapoly)), rhoSample),'r');
        title('fitted dTa/drho')
    elseif AlphaType==2
        figure;
        plot(rhoSample,PaSample,'b');
        hold on;
        plot(rhoSample,polyval(fliplr(Papoly), rhoSample),'r');
        legend('original Pa', 'fitted Pa'); 
        title('Pa')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Papoly)), rhoSample),'r');
        title('fitted dPa/drho')
    end
end


if BeamType~=3
    figure;
    plot(rhoSample,nbSample,'b');
    hold on;
    plot(rhoSample,polyval(fliplr(nbpoly), rhoSample),'r');
    legend('original nb', 'fitted nb'); 
    title('nb')
    figure;
    plot(rhoSample,polyval(polyder(fliplr(nbpoly)), rhoSample),'r');
    title('fitted dnb/drho')
    if BeamType==1
        figure;
        plot(rhoSample,TbSample,'b');
        hold on;
        plot(rhoSample,polyval(fliplr(Tbpoly), rhoSample),'r');
        legend('original Tb', 'fitted Tb'); 
        title('Tb')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Tbpoly)), rhoSample),'r');
        title('fitted dTb/drho')
    elseif BeamType==2
        figure;
        plot(rhoSample,PbSample,'b');
        hold on;
        plot(rhoSample,polyval(fliplr(Pbpoly), rhoSample),'r');
        legend('original Pb', 'fitted Pb'); 
        title('Pb')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Pbpoly)), rhoSample),'r');
        title('fitted dPb/drho')
    end
end