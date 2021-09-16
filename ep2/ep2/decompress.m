function [] = decompress(compressedImg, method, k, h)
  filename = "decompressed.png";
  n = 200;
  
  # definições de x barra e y barra
  x = 0;
  y = 0;
  
  # PARTE 1 - O ZOOLOGICO
  if (compressedImg == -1)
    ## constroi a imagem a partir da função dada
    for (i = 1:n)
      for (j = 1:n)
        imagem(i,j,:) = f(x + i*h,y + j*h);
      endfor
    endfor
    imwrite(imagem,"decompressAntes.png");
  # PARTE 2 - A SELVA
  else
    ## carrega a imagem
    imagem = imread(compressedImg);
    imagem = im2double(imagem);
  endif
  
  ## método que será aplicado para construir a nova imagem
  if (method == 1) # bilinear
    imagem = bilinear(imagem, k, h);
  elseif (method == 2) # bicubico  
    imagem = bicubica(imagem, k, h);
  endif
  
  imwrite(imagem,filename,'Compression','none');
endfunction

function num = f(x, y)
  num = [sin(x)+cos(y),sin(x)+cos(y),sin(x)+cos(y)];
endfunction

function novaImagem = bilinear(matrizImagem, k, h)
  n = length(matrizImagem);
  # seta o novo h para que ele seja n+k partes iguais do intervalo procurado
  # que era n*h
  h_linha = n*h/(n +(n-1)*k)
  
  proporcao = ((n+(n-1)*k)/n);
  
  H = [1,0,0,0;
       1,0,h,0;
       1,h,0,0;
       1,h,h,h^2];
  
  # constroi a nova imagem
  for cor = 1:3
    for ii = 1:(n + (n-1)*k)
      for jj = 1:(n +(n-1)*k)
        # esse é a posição do vertice do quadrado Qij que vamos verificar
        # cujo (x,y) esta dentro
        i = floor(ii/proporcao)+1;
        j = floor(jj/proporcao)+1;
        
        # definimos x e y
        x = jj*h_linha;
        y = ii*h_linha;
        
        # define o vetor [f(xi,yi), f(xi,yi+1), f(xi+1,yi), f(xi+1,yi+1)]
        if (i >= n && j >= n)
          F = [matrizImagem(n-1,n-1,cor);
               matrizImagem(n-1,n,cor);
               matrizImagem(n,n-1,cor);
               matrizImagem(n,n,cor)];
        elseif (i >= n)
          F = [matrizImagem(n-1,j,cor);
               matrizImagem(n-1,j+1,cor);
               matrizImagem(n,j,cor);
               matrizImagem(n,j+1,cor)];
        elseif (j >= n)
          F = [matrizImagem(i,n-1,cor);
               matrizImagem(i,n,cor);
               matrizImagem(i+1,n-1,cor);
               matrizImagem(i+1,n,cor)];
        else
          F = [matrizImagem(i,j,cor);
               matrizImagem(i,j+1,cor);
               matrizImagem(i+1,j,cor);
               matrizImagem(i+1,j+1,cor)];
        endif
               
          A = H\F;
          novaImagem(ii,jj,cor) = A'(1) + A'(2)*(x - j*h) + A'(3)*(y - i*h) + A'(4)*(x - j*h)*(y - i*h);
      endfor
    endfor
  endfor
endfunction

function novaImagem = bicubica(matrizImagem, k, h)
  n = length(matrizImagem);
  # seta o novo h para que ele seja n+k partes iguais do intervalo procurado
  # que era n*h
  h_linha = n*h/(n +(n-1)*k)
  
  proporcao = ((n+(n-1)*k)/n);
  
  B = [1,0,0,0;
       1,h,h^2,h^3;
       0,1,0,0;
       0,1,2*h,3*h^2];
  
  # constroi a nova imagem
  for cor = 1:3
    for ii = 1:(n + (n-1)*k)
      for jj = 1:(n +(n-1)*k)
        # esse é a posição do vertice do quadrado Qij que vamos verificar
        # cujo (x,y) esta dentro
        i = floor(ii/proporcao)+1;
        j = floor(jj/proporcao)+1;
        
        # definimos x e y
        x = jj*h_linha;
        y = ii*h_linha;
        
        # matriz dos pontos que formam Qij
        if (i >= n && j >= n)
          F = [matrizImagem(n-1,n-1,cor),matrizImagem(n-1,n,cor);
               matrizImagem(n,n-1,cor),  matrizImagem(n,n,cor)];
        elseif (i >= n)
          F = [matrizImagem(n-1,j,cor),  matrizImagem(n-1,j+1,cor);
               matrizImagem(n,j,cor),    matrizImagem(n,j+1,cor)];
        elseif (j >= n)
          F = [matrizImagem(i,n-1,cor),  matrizImagem(i,n,cor);
               matrizImagem(i+1,n-1,cor),matrizImagem(i+1,n,cor)];
        else
          F = [matrizImagem(i,j,cor),    matrizImagem(i,j+1,cor);
               matrizImagem(i+1,j,cor),  matrizImagem(i+1,j+1,cor)];
        endif
        
        # derivadas primeiras em relação a x e y
        if (j > n)
          dxF = achadxF(matrizImagem, i, n, n, h, cor);
        else
          dxF = achadxF(matrizImagem, i, j, n, h, cor);
        endif
        
        if (i > n)
          dyF = achadyF(matrizImagem, n, j, n, h, cor);
        else
          dyF = achadyF(matrizImagem, i, j, n, h, cor);
        endif
        
        ## calcula as derivadas segundas em relação a y
        for numI = 1:2
          for numJ = 1:2
            dxdyF(numI,numJ) = (dyF(numI+2,numJ) - dyF(numI,numJ))/(2*h);
          endfor
        endfor

        #H  = [F,dyF(2:3,:);dxF(2:3,:),dxdyF]
        
        A = B\[F,dyF(2:3,:);dxF(2:3,:),dxdyF]/(B');
        
        xi = j*h;
        yi = i*h;
        
        novaImagem(ii,jj,cor) = [1,(x-xi),(x-xi)^2,(x-xi)^3]*A*[1;(y-yi);(y-yi)^2;(y-yi)^3];
      endfor
    endfor
  endfor
endfunction

function dxF = achadxF(matrizImagem, i, j, n, h, cor)
  if (i == 1)
    for numI = 1:4
      for numJ = 1:2
        ## calcula as derivadas primeiras em relação a x
        dxF(numI,numJ) = (matrizImagem(2,j,cor) - matrizImagem(1,j,cor))/(2*h);
      endfor
    endfor
  elseif (i >= n)
    for numI = 1:4
      for numJ = 1:2
        ## calcula as derivadas primeiras em relação a x
        dxF(numI,numJ) = (matrizImagem(n,j,cor) - matrizImagem(n-1,j,cor))/(2*h);
      endfor
    endfor
  else
    for numI = 1:4
      for numJ = 1:2
        ## calcula as derivadas primeiras em relação a x
        dxF(numI,numJ) = (matrizImagem(i+1,j,cor) - matrizImagem(i-1,j,cor))/(2*h);
      endfor
    endfor
  endif
endfunction

function dyF = achadyF(matrizImagem, i, j, n, h, cor)
  if (j == 1)
    for numI = 1:4
      for numJ = 1:2
        ## calcula as derivadas primeiras em ralação a y
        dyF(numI,numJ) = (matrizImagem(i,2,cor) - matrizImagem(i, 1, cor))/(2*h);
      endfor
    endfor
  elseif (j >= n)
    for numI = 1:4
      for numJ = 1:2
        ## calcula as derivadas primeiras em ralação a y
        dyF(numI,numJ) = (matrizImagem(i, n, cor) - matrizImagem(i, n-1, cor))/(2*h);
      endfor
    endfor
  else
    for numI = 1:4
      for numJ = 1:2
        ## calcula as derivadas primeiras em ralação a y
        dyF(numI,numJ) = (matrizImagem(i,j+1,cor) - matrizImagem(i,j-1,cor))/(2*h);
      endfor
    endfor
  endif
endfunction