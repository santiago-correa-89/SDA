function file = defineFile(t, path)

    if strcmp(path, 'data/Tjaereborg')
        
        if t == 100
            file= sprintf('AirfoilData/tjaereCil.dat');
        else
            file = sprintf('AirfoilData/tjaere%02d_ds.dat', t);
        end
    
    elseif strcmp(path, 'Data/NREL2.5MW_116')
        
        if t == 100
            file = sprintf('AirfoilData/airfoil000_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 97
            file = sprintf('AirfoilData/airfoil001_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 94
            file = sprintf('AirfoilData/airfoil002_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 91
            file = sprintf('AirfoilData/airfoil003_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 88
            file = sprintf('AirfoilData/airfoil004_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 85
            file = sprintf('AirfoilData/airfoil005_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 82
            file = sprintf('AirfoilData/airfoil006_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 79
            file = sprintf('AirfoilData/airfoil007_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 76
            file = sprintf('AirfoilData/airfoil008_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 73
            file = sprintf('AirfoilData/airfoil009_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 70
            file = sprintf('AirfoilData/airfoil010_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 67
            file = sprintf('AirfoilData/airfoil011_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 64
            file = sprintf('AirfoilData/airfoil012_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 61
            file = sprintf('AirfoilData/airfoil013_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 58     
            file = sprintf('AirfoilData/airfoil014_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 55
            file = sprintf('AirfoilData/airfoil015_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 52     
            file = sprintf('AirfoilData/airfoil016_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 49     
            file = sprintf('AirfoilData/airfoil017_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 46     
            file = sprintf('AirfoilData/airfoil018_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 43     
            file = sprintf('AirfoilData/airfoil019_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 40     
            file = sprintf('AirfoilData/airfoil020_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 37     
            file = sprintf('AirfoilData/airfoil021_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 34     
            file = sprintf('AirfoilData/airfoil022_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 31     
            file = sprintf('AirfoilData/airfoil023_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 28     
            file = sprintf('AirfoilData/airfoil024_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 25     
            file = sprintf('AirfoilData/airfoil025_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 22    
            file = sprintf('AirfoilData/airfoil026_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 19     
            file = sprintf('AirfoilData/airfoil027_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 16     
            file = sprintf('AirfoilData/airfoil028_XXXXXXX_%03d.wtgAirfoil', t);
        elseif t == 13     
            file = sprintf('AirfoilData/airfoil029_XXXXXXX_%03d.wtgAirfoil', t);
        end

    else
        
        error('Ruta no reconocida.');
    
    end

end