function mu = water(x,y,z)
  r2 = ((x-0.05).*(x-0.05) + (y-0.05).*(y-0.05) + z.*z);
  mu = ones(size(r2)).*4;

  for i=1:length(r2(:))
    if ( r2(i) < 0.36)
      mu(i) = 80;
    elseif (r2(i) > 0.49)
      mu(i) = 4;
    else
      mu(i) = 80 - (80.0 - 4.0)/(0.49 - 0.36)*(r2(i) - 0.36);
    end
  end % for
end 