function a = mkchoice(p)

a = max(find([-eps cumsum(p)] < rand));

return