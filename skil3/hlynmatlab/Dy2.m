function y=Dy2(y1,y2)
    gamma =0.005; delta=0.2;
    y=gamma.*y1.*y2-delta.*y2;
end
