function dh = load_dh_parameters
    dh.d = [0.089159 0 0 0.10915 0.09465 0.0823];
    dh.m = [3.7000 8.3930 2.2750 1.2190 1.2190 0.1879];
    dh.alpha = [pi/2 0 0 pi/2 -pi/2 0];
    dh.a = [0 -0.42500 -0.39225 0 0 0];
    cm1 = [0.0 -0.02561 0.00193];
    cm2 = [0.2125 0.0 0.11336];
    cm3 = [0.15 0.0 0.0265];
    cm4 = [0.0 -0.0018 0.01634];
    cm5 = [0.0 0.0018 0.01634];
    cm6 = [0.0 0.0 -0.001159];

    dh.cm = [cm1' cm2' cm3' cm4' cm5' cm6'];
    dh.dof = 6;
    dh.issym = false;

    i1 = [0.010267 0.010267 0.00666];
    i2 = [0.2269 0.2269 0.0151];
    i3 = [0.0312168 0.0312168 0.004095];
    i4 = [0.002559898976 0.002559898976 0.0021942];
    i5 = [0.002559898976 0.002559898976 0.0021942];
    i6 = [8.46958911216e-5 8.46958911216e-5 0.0001321171875];

    dh.I = zeros(3, 3, 6);
    dh.I(:,:,1) = [i1(1) 0 0; 0 i1(2) 0; 0 0 i1(3)];
    dh.I(:,:,2) = [i2(1) 0 0; 0 i2(2) 0; 0 0 i2(3)];
    dh.I(:,:,3) = [i3(1) 0 0; 0 i3(2) 0; 0 0 i3(3)];
    dh.I(:,:,4) = [i4(1) 0 0; 0 i4(2) 0; 0 0 i4(3)];
    dh.I(:,:,5) = [i5(1) 0 0; 0 i5(2) 0; 0 0 i5(3)];
    dh.I(:,:,6) = [i6(1) 0 0; 0 i6(2) 0; 0 0 i6(3)];

    gravity = [0 0 9.81]';
end