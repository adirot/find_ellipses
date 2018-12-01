function isInsideEllipse = isInsideEllipse(x,y,ellipseParam)
    X0 = ellipseParam.X0_in;
    Y0 = ellipseParam.Y0_in;
    a = ellipseParam.a;
    b = ellipseParam.b;
    phi = ellipseParam.phi;
    isInsideEllipse=logical(((((((cos(phi).*(x-X0))+(sin(phi).*(y-Y0)))/a).^2)+((((sin(phi).*(x-X0))-(cos(phi).*(y-Y0)))/b).^2)) <= 1));
end