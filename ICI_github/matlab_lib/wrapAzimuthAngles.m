function phi = wrapAzimuthAngles(phi)

    phi =  (mod(phi + 180,360) - 180);

end