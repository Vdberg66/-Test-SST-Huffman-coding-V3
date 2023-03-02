%              f_Statistics_Set_Shaping_Theory_Huffman_l
%                        Contact info
%                Christian.Schmidt55u@gmail.com
%                    adrain.vdberg66@yahoo.com 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Definition
%
% Definition: Given a sequence x we call Cs(x) the coded sequence in which 
% the symbols are replaced by a uniquely decodable code. 
%
% Definition: given a sequence x of random variables of length N, we call 
% the coding limit Lc(x) the function defined by the following example:
% 
% given a sequence:
% 123445
%  
% We calculate the frequencies of the symbols present in the sequence.
%
%  symbol 1  frequencies 1/6
%  symbol 2  frequencies 1/6
%  symbol 3  frequencies 1/6
%  symbol 4  frequencies 2/6
%  symbol 5  frequencies 1/6
%
%  Lc(x)=-log2(1/6)-log2(1/6)-log2(1/6)-log2(2/6)-log2(2/6)-log2(1/6)=13.5 bit
%
% 13.5 bit represents the minimum length of the coded sequence 
% in which the symbols have been replaced with uniquely decodable code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The program performs the following operations:
% 1) generates a random sequence x with uniform distribution
% 2) calculate the frequencies of the symbols present in the sequence x 
% 3) use this information to calculate the coding limit Lc(x)
% 4) apply the transform f(x)
% 5) code the transformed sequence
% 6) compares the coding limit Lc(x) of the generated sequence 
%    with the length of the encoded transformated sequence cs(f(x)) 
% 7) repeats all these steps a number of times defined by the parameter history
% 8) display the average values obtained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Important
%
% If you change the length of the sequence and the number of the generated 
% symbols, you have to be careful that the Huffman encoding approximates the 
% coding limit of the sequence by about one bit. if you take too 
% long sequences the Huffman algorithm becomes very inefficient therefore, 
% it cannot detect the advantage obtained by the transformation.
% As a general rule, if you take ns symbols the length of the sequence 
% must be about 2*ns, in this case the Huffman encoding approximates 
% the coding limit of the sequence by about one bit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Results for different settings
%
%  The fSST2 and invfSST2 function only works when the number of symbols ns 
%  is greater or equal than 20 and less or equal than 500 and the length of 
%  the sequence N is greater or equal than 40 and less or equal than 1000. 
%  So, it is recommended to generate random sequences with a number of 
%  symbols and length greater than Ns=20 and N=40.
%
%  Ps=probability with which the transformed sequence f(x) can be encoded 
%  using a uniquely decodable code (Huffman coding) with lower bit number 
%  than coding limit Lc(x) of the initial sequence x.
%
%  ns= number of symbols
%  N=length of the sequence
%
%  ns   N    Ps 
%  40   80   69%
%  50  100   72%
%  60  120   80%
% 500 1000   88%
%
%  As you can see, increasing the symbol number increases the probability Ps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
history=100;
ns=32;
len=64;
cs=0;
totcodel=0;
totlc=0;
tottlc=0;
itnent=0;
lc=0;
tlc=0;
itlc=0;

for i=1:history

 % Generation of the sequence with a uniform distribution

 symbols=1:ns;
 prob(1,1:ns)=1/ns;
 seq=randsrc(1,len,[symbols; prob]);

 % coding limit Lc(x)

 lc=0;

 for i2=1:len

  sy=seq(1,i2);
  fs=nnz(seq==sy)/len;
  lc=lc-log2(fs);

 end

 % Start trasformation

 mcodel=10000;

 nseq=fSSTt(seq);

 % The new sequence is long nlen=len+1
 
 nlen=len+1;

 % coding limit of the transformed sequence of length nlen

  tlc=0;

  for i2=1:nlen

   sy=nseq(1,i2);
   fs=nnz(nseq==sy)/nlen;
   tlc=tlc-log2(fs);

  end 

 % Having transformed the sequence, we have to redefine the length of the 
 % vectors that are used in the encoding

 index=0;

 for i2=1:ns

  fs=nnz(nseq==i2)/nlen;

  if fs > 0

   index=index+1;    

  end

 end 

 c=zeros(index,1);
 vs=zeros(index,1);
 index=0;

 % We calculate the frequencies of the symbols in the transformed sequence

 for i2=1:ns

  fs=nnz(nseq==i2)/nlen;

  if fs > 0

   index=index+1;    
   c(index)=nnz(nseq==i2)/nlen;  
   vs(index)=i2;

  end

 end 

 % We code the transformed sequence

 counts=c;

 tdict=huffmandict(vs,counts);
 tcode=huffmanenco(nseq,tdict);

 bcode=de2bi(tcode);
 tcodel=numel(bcode);

 % If the length of the encoded message of the transformed sequence is less
 % than coding limit Lc(x) of the original sequence x, we increase the counter
 % cs by one

 if tcodel < lc

  cs=cs+1;

 end

 % We apply the inverse transform and we obtain the initial sequence

 iseq=invfSSTt(nseq);

 % We check that the obtained sequence is equal to the initial sequence

 flag=isequal(seq,iseq);

 if flag == false

   fprintf('Error, sequence not equal to the initial sequence\n');

 end

 totcodel=totcodel+tcodel;
 totlc=totlc+lc;
 tottlc=tottlc+tlc;

end

% We calculate the average of the coding limit Lc(x) of the generated sequences
% x,the average of the coding limit Lc(f(x)) of the transformed sequences f(x) 
% and the average of the length of the encoded transformed sequence Cs(f(x))
  
 medlc=totlc/history;
 medcodel=totcodel/history;
 medtlc=tottlc/history;

% We calculate the percentage of sequences where the length of the encoded
% transformed sequence Cs(f(x)) is less than the coding limit Lc(x) of the 
% generated sequence x

 pcs=(cs/history)*100;

% We display the average values obtained

 fprintf('The average of the coding limit Lc(x) of the generated sequences x\n');
 medlc

 fprintf('The average of the length of the encoded transformed sequence Cs(f(x)) \n');
 medcodel

 fprintf('The average of the coding limit Lc(f(x)) of the transformed sequences f(x)\n');
 medtlc

 fprintf('Number of sequences where the length of the encoded transformed sequence Cs(f(x)) is less than the coding limit Lc(x) of the generated sequence x\n');
 cs

 fprintf('There is a percentage of %2.0f%% that length of the encoded transformed sequence Cs(f(x)) is less than the coding limit Lc(x) of the generated sequence x\n',pcs);

