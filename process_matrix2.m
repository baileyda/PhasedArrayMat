function[matrixmag,matrixdB,matrixdBnorm] = process_matrix2(matrix)
matrixmag=abs(matrix);
matrixmax=max(max(matrixmag));
matrixdB=10*log10(matrixmag+eps);
matrixdBnorm=matrixdB - 10*log10(matrixmax);
