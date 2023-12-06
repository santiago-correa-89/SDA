function datos = importAirfoilData(path, file)

try
    if strcmp(path, 'Data/NREL2.5MW_116')
        datos = importdata(fullfile(path, file));
    elseif strcmp(path, 'data/Tjaereborg')
        datos = importdata(fullfile(path, file));
    end
catch
    disp('No se pudo cargar el archivo.');
end

end
