function y=Dy1(y1,y2)
    alpha=0.5; beta=0.01;
    y=alpha.*y1-beta.*y1.*y2;
end
    