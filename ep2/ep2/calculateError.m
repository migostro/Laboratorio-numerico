function err = calculateError(originalImg, decompressedImg)
  orig = imread(originalImg);
  orig = im2double(orig);
  
  dec = imread(decompressedImg);
  dec = im2double(dec);
  
  errR = (norm(orig(:,:,1)) - norm(dec(:,:,1)))^2/((norm(orig(:,:,1)))^2);
  errG = (norm(orig(:,:,2)) - norm(dec(:,:,2)))^2/((norm(orig(:,:,2)))^2);
  errB = (norm(orig(:,:,3)) - norm(dec(:,:,3)))^2/((norm(orig(:,:,3)))^2);
  err = (errR+errG+errB)/3;
endfunction
