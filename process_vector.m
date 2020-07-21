function[vectormag,vectordB,vectordBnorm] = process_vector(vector)
vectormag=abs(vector);
vectordB=20*log10(vectormag+eps);
vectordBnorm=20*log10((vectormag+eps)/max(vectormag));
