fields := [2, 3, 4, 5, 7, 8, 9, 11];

// Express a form of bidegree (d,d) in Segre coordinates (TU, TV, SU, SV).
function SegreForm(form, d, x)
    if form eq 0 then
        return 0*x[1];
    end if;
    result := 0*x[1];
    monomials := Monomials(form);
    coefficients := Coefficients(form);
    for n in [1..#monomials] do
        i := Degree(monomials[n], 1);
        j := Degree(monomials[n], 3);
        r := Min(i, j);
        result +:= coefficients[n] * x[1]^r * x[2]^(i-r) * x[3]^(j-r) * x[4]^(d-Max(i,j));
    end for;
    return result;
end function;

function Example(q)
    assert q in {2, 3, 4, 5, 7, 8, 9, 11};
    k := GF(q);
    base<T,S,U,V> := ProductProjectiveSpace(k, [1, 1]);
    weighted<x0,x1,x2,x3,w> := ProjectiveSpace(k, [1, 1, 1, 1, 2]);


    if q eq 2 then
        H := (T^2 + T*S + S^2)*(U^2 + U*V + V^2);
        F := H^2 + T*S*(T + S)*(S*U^3*V + T*U*V^3);

    elif q eq 3 then
        H := 0;
        A := T*S*(T^2 + S^2);
        B := 2*(T + S)*(T^3 - T^2*S + S^3);
        Q1 := U^2 + V^2;
        Q2 := U^2 + U*V - V^2;
        F := B*Q1^2 - A*Q2^2;

    elif q eq 4 then
        a := k.1;
        A := a*(S*U + T*V) + S*V;
        B := (1 + a)*(T*U + T*V) + a*S*U;
        H := A^2 + A*B + a*B^2;
        F := (T^4 + a*(S^4 + S*T^3))*U^4
             + (1 + a)*(T^4 + S^4 + S*T^3)*(U^3*V + V^4);

    elif q eq 5 then
        H := 0;
        A := 2*S*T*(T - S)*(T - 2*S);
        B := 2*(T - 3*S)*(T - 4*S)*(T^2 + 3*S^2);
        Q1 := U^2 + 2*V^2;
        Q2 := U^2 + U*V + V^2;
        F := B*Q1^2 - A*Q2^2;

    elif q eq 7 then
        H := 0;
        A := T*(T - S)*(T - 2*S)*(T - 4*S);
        B := 3*S*(T - 3*S)*(T - 5*S)*(T - 6*S);
        Q1 := U^2 + V^2;
        Q2 := 2*U^2 - U*V - 2*V^2;
        F := B*Q1^2 - A*Q2^2;

    elif q eq 8 then
        a := k.1;
        A := (1 + a + a^2)*(T*U + S*V) + (1 + a^2)*T*V;
        B := (1 + a)*T*U + (a + a^2)*S*U + a*S*V;
        H := A^2 + A*B + B^2;
        F := (T^4 + T*S^3 + S^4)*U^4
             + (T^4 + T*S^3 + (1 + a)*S^4)*U*V^3
             + ((1 + a^2)*T^4 + (1 + a + a^2)*T*S^3 + S^4)*V^4;

    elif q eq 9 then
        i := Sqrt(k!-1);
        d := i - 1;
        Q := U^2 + d*V^2;
        Z := U^2 - d*V^2;
        H := 0;
        F := (1 + i)*(S^4 + i*T^4)*Q^2
             - T*S*Z*(T*S*Z + i*U*V*(S^2 - d*T^2));

    elif q eq 11 then
        H := 0;
        A := S*T*(T - S)*(T - 2*S);
        B := (T - 3*S)*(T - 4*S)*(T - 6*S)*(T - 8*S);
        Q1 := U^2 - 2*U*V - V^2;
        Q2 := U^2 + 2*U*V - V^2;
        F := 6*(B*Q1^2 - A*Q2^2);
    end if;

    x := [x0, x1, x2, x3];
    segre := SegreEmbedding(base);
    base_equation := Evaluate(DefiningPolynomial(segre), x);
    double_cover := Scheme(weighted, [base_equation, w^2 + SegreForm(H, 2, x)*w - SegreForm(F, 4, x)]);
    segre_embedded := Image(map<double_cover -> ProjectiveSpace(k, 10) | Sections(LinearSystem(weighted, 2))>);
    return double_cover, segre_embedded;
end function;

for q in fields do
    double_cover, segre_embedded := Example(q);
    assert IsNonsingular(double_cover);
    assert #RationalPoints(segre_embedded) eq 0;
    printf "F_%o: smooth and pointless\n", q;
end for;
