rhoSample = linspace(0,1,512);
neSample = polyval(fliplr(nepoly), rhoSample);
TeSample = polyval(fliplr(Tepoly), rhoSample);

figure;
plot(rhoSample,neSample,'r');
title('ne')
figure;
plot(rhoSample,polyval(polyder(fliplr(nepoly)), rhoSample),'r');
title('dne/drho')

figure;
plot(rhoSample,TeSample,'r');
title('Te')
figure;
plot(rhoSample,polyval(polyder(fliplr(Tepoly)), rhoSample),'r');
title('dTe/drho')

if IonType~=3

    niSample = polyval(fliplr(nipoly), rhoSample);
    figure;
    plot(rhoSample,niSample,'r');
    title('ni')
    figure;
    plot(rhoSample,polyval(polyder(fliplr(nipoly)), rhoSample),'r');
    title('dni/drho')

    if IonType==1

        TiSample = polyval(fliplr(Tipoly), rhoSample);
        figure;
        plot(rhoSample,TiSample,'r');
        title('Ti')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Tipoly)), rhoSample),'r');
        title('dTi/drho')

    elseif IonType==2

        PiSample = polyval(fliplr(Pipoly), rhoSample);
        figure;
        plot(rhoSample,PiSample,'r');
        title('Pi')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Pipoly)), rhoSample),'r');
        title('dPi/drho')

    end

end

if AlphaType~=3

    naSample = polyval(fliplr(napoly), rhoSample);
    figure;
    plot(rhoSample,naSample,'r');
    title('na')
    figure;
    plot(rhoSample,polyval(polyder(fliplr(napoly)), rhoSample),'r');
    title('dna/drho')

    if AlphaType==1

        TaSample = polyval(fliplr(Tapoly), rhoSample);
        figure;
        plot(rhoSample,TaSample,'r');
        title('Ta')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Tapoly)), rhoSample),'r');
        title('dTa/drho')

    elseif AlphaType==2

        PaSample = polyval(fliplr(Papoly), rhoSample);
        figure;
        plot(rhoSample,PaSample,'r');
        title('Pa')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Papoly)), rhoSample),'r');
        title('dPa/drho')

    end
    
end

if BeamType~=3

    nbSample = polyval(fliplr(nbpoly), rhoSample);
    figure;
    plot(rhoSample,nbSample,'r');
    title('nb')
    figure;
    plot(rhoSample,polyval(polyder(fliplr(nbpoly)), rhoSample),'r');
    title('dnb/drho')

    if BeamType==1

        TbSample = polyval(fliplr(Tbpoly), rhoSample);
        figure;
        plot(rhoSample,TbSample,'r');
        title('Tb')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Tbpoly)), rhoSample),'r');
        title('dTb/drho')

    elseif BeamType==2

        PbSample = polyval(fliplr(Pbpoly), rhoSample);
        figure;
        plot(rhoSample,PbSample,'r');
        title('Pb')
        figure;
        plot(rhoSample,polyval(polyder(fliplr(Pbpoly)), rhoSample),'r');
        title('dPb/drho')

    end
    
end