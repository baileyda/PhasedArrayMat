function[matrixmag,matrixdB,matrixdBnorm] = process_matrix(matrix)
matrixmag=abs(matrix);
matrixmax=max(max(matrixmag));
matrixdB=20*log10(matrixmag+eps);
matrixdBnorm=matrixdB - 20*log10(matrixmax);
