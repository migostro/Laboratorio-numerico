function [] = compress(originalImg, k)
  nome_arquivo = "compressed.png";
  imagem = imread(originalImg,"png");
  n = length(imagem);
  for cor = 1:3
    i = 1;
    for ii = 1:(n)
      j = 1;
      for jj = 1:(n)
        if (mod(ii,k+1) == 0 && mod(jj,k+1) == 0)
          novaImagem(i,j,cor) = imagem(ii,jj,cor);
          j = j+1;
        endif
      endfor
      if (mod(ii,k+1) == 0)
        i = i+1;
      endif
    endfor
  endfor
  imwrite(novaImagem, nome_arquivo,"png");
endfunction