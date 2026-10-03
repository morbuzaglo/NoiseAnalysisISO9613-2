function D_est = estimate_Dsource_at_angle(anglesDeg,D_source,thetaDeg)

    anglesDeg = mod(anglesDeg(:),360);
    thetaDeg = mod(thetaDeg,360);

    [anglesDeg,idx] = sort(anglesDeg);
    D_source = D_source(idx,:);

    anglesExt = [anglesDeg(end)-360; anglesDeg; anglesDeg(1)+360];
    D_ext = [D_source(end,:); D_source; D_source(1,:)];

    D_est = interp1(anglesExt,D_ext,thetaDeg,'linear');

end