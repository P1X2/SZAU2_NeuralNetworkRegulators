run("uczenie.m")

if farx(length(farx)-1) < bestARX
    bestARX = farx(length(farx)-1);
end

if foe(length(foe)-1) < bestOE
    bestOE = foe(length(foe)-1);
    run("model.m")
end