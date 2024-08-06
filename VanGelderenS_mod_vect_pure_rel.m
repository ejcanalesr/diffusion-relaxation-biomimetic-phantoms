function Svg = VanGelderenS_mod_vect_pure_rel(x, Data)
    r     = x(1);
    beta  = x(2);
    % Data = [TE, T2csf, rho_j])
    rho_j = Data(end); Data(end) = [];
    T2csf = Data(end); Data(end) = [];
    TE    = Data;
    Svg = beta * exp(-TE./T2csf) .* exp(-2*rho_j*TE./r);
end