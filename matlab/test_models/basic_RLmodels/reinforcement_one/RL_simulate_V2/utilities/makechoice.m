function a = makechoice(p)

% choose an option based on the choice current probabiltiies
a = max(find([-eps cumsum(p)] < rand));

end 