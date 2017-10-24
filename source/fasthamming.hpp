#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <string.h>
static char* reverseBytes(void *ptr, char nBytes);
static void printBits(void *ptr, char nBytes);
void printBitsEndian(void *ptr, char nBytes, char endian); 
 //Encoding is 
 //0000 N  0
 //0001 A  1
 //0010 C  2
 //0100 G  4
 //1000 T  8
 
//fastest is to use 64 bit version of popcount
//32 bit LUT is second fastest - use when popcount not available
//function pointer to choose is too slow - use ifdefs
//inlining is fastest

//Panel case
//case when there are no N's in one sequence (panel sequence)
//simplifies calculations because there is no NN and all NX comes from other sequence
//bitcount (a^b) gives 2XY+NX;
//multiply by 2 and add number of Ns to get correct distance of  4*distance=4*XY+ 3*NX

//General case where there are uncertain bases in both sequences is more difficult
//need to keep track of Ns in separate and correct the distance
//first estimate is bitcount(a^b);
//then correct for NX and NN
//NN=bitcount(na & nb)
//correction factor is 16*NN+15*(2*NN-na.Ncount-nb.Ncount) added to 8*distance estimate to give 16*distance

using namespace std;

template <class T> class Nvector{
	//bitarray for storing N positions
	//T is an integer class - usually 64 bit for pop sequence
	public:
 	vector <unsigned char> sequence;
  size_t nBytes;
  size_t Ncount;
  T* Tseq; //use this for comparison
  Nvector(){
		 sequence.resize(0);
		 Ncount=0;
		 Tseq=(T*) &(sequence[0]);
		 nBytes=0;			
		}	
  Nvector(size_t nBits){
			const int Tbits=sizeof(T)*8;
		 nBytes=(nBits%Tbits) ? (nBits/Tbits+1)*sizeof(T): (nBits/Tbits)*sizeof(T);
		 //pad nBytes;
		 sequence.resize(nBytes);
		 Ncount=0;
		 memset(&(sequence[0]),0,nBytes);
		 Tseq=(T*) &(sequence[0]);
	 }
	 Nvector(const Nvector &A){
  	nBytes=A.nBytes;
	 	Ncount=A.Ncount;
	  sequence=A.sequence;
	  Tseq=(T*)&(sequence[0]);
	 }
	 Nvector & operator=( const Nvector &A ) {
   nBytes=A.nBytes;
	 	Ncount=A.Ncount;
	  sequence=A.sequence;
	  Tseq=(T*)&(sequence[0]);
   return *this;
		}
	 unsigned char getBit(size_t index) {
   return (sequence[index/8] >> 7-(index & 0x7)) & 0x1;
  }
 void setBit(size_t index) {
		if(index/8 >= nBytes){
			nBytes*=2;
			sequence.resize(nBytes);
			Tseq=(T*) &(sequence[0]);
		}
  sequence[index/8] = sequence[index/8] |  1 << 7-(index & 0x7);
  Ncount++;
 }
 void reset(){
		memset(&(sequence[0]),0,nBytes);
		Ncount=0;
	}	
};
template <class T> class acgt{
	//Encoding is 
 //0000 N  0
 //0001 A  1
 //0010 C  2
 //0100 G  4
 //1000 T  8
	//property being that hamming bitwise distance between any two bases ACGT is 2 an is 1 to N
	public:
		std::vector <unsigned char> sequence; //padded to sizeofT bytes
  uint32_t seqLength;
  unsigned char *cseq; //aliases to sequence
  T *Tseq; //aliases to sequence
  unsigned int Ncount=0;
	 acgt(const char *inputSequence,int _len){
			if(! _len) _len=strlen(inputSequence);
		 seqLength=_len;
		 const int nBytes=(seqLength%(2*sizeof(T)))? (seqLength/(2*sizeof(T))+1)*sizeof(T) : seqLength/(2*sizeof(T))*sizeof(T); //bytes necessary to hold sequence
   sequence.resize(nBytes);
   memset(&(sequence[0]),0,sequence.size());
		 cseq=(unsigned char*) &(sequence[0]);
		 Tseq=(T*)  &(sequence[0]);
		 seqLength=seqLength;
	  for (int i=0;i<seqLength;i++){
		  uint8_t k=(i%2)? 1 : 1<<4;
		  if(inputSequence[i]=='A'){
		 	cseq[i/2]|=k;
		 	continue;
		 }
		 k=k<<1;	
		 if(inputSequence[i]=='C'){
		  cseq[i/2]|=k;
    continue;
		 }
		 k=k<<1;	
		 if(inputSequence[i]=='G'){
		 	cseq[i/2]|=k;
    continue;
		 }
		 k=k<<1;
		 if(inputSequence[i]=='T'){
		 	cseq[i/2]|=k;
		 	continue;
		 }
		 else Ncount++;
		}
	}
	acgt(){
		seqLength=0;
  sequence.resize(0);
		cseq=(unsigned char*) &(sequence[0]);
		Tseq=(T* ) &(sequence[0]);
		Ncount=0;
	}	
	acgt(int len){
		seqLength=len;
  const int nElements=(seqLength%(2*sizeof(T)))? (seqLength/(2*sizeof(T))+1) : seqLength/(2*sizeof(T)); //bytes necessary to hold sequence
  const int nBytes=nElements*sizeof(T);
  sequence.resize(nBytes);
		memset(&(sequence[0]),0,sequence.size());
		cseq=(unsigned char*) &(sequence[0]);
		Tseq=(T* ) &(sequence[0]);
		Ncount=0;
	}
		
	acgt(const acgt &a){
		sequence=a.sequence;
		cseq=(unsigned char*) &(sequence[0]);
		Tseq=(T*) &(sequence[0]);
		seqLength=a.seqLength;
		Ncount=a.Ncount;
	}
	acgt & operator=( const acgt&a ) {
  sequence=a.sequence;
		cseq=(unsigned char*) &(sequence[0]);
		Tseq=(T*) &(sequence[0]);
		seqLength=a.seqLength;
		Ncount=a.Ncount;
  return *this;
 }
	acgt subsequence(int start, int len){
		acgt sub(len);
		sub.Ncount=0;
		if(!Ncount){
   if(start % 2){
		 	//have to copy byte by byte and shift left by 4;
		 	int k=start/2;
		  for(int i=0;i<len/2;i++){
		 		sub.cseq[i]=(cseq[k] << 4) | (cseq[k+1] >>4);
		 		k++;
			 }
			 //don't forget extra half byte	
		 }
		 else{	
    memcpy(sub.cseq,cseq+start/2,len/2);
		 }
		 //don't forget to copy odd byte
	 	if(len%2){
	 		sub.cseq[len/2]=cseq[start/2+len/2] & 0xF0;
	 	}
		}
		else{
			 //count at the same time
		  if(start % 2){
		 	//have to copy byte by byte and shift left by 4;
		 	int k=start/2;
		  for(int i=0;i<len/2;i++){
					unsigned char upper=cseq[k] << 4;
					unsigned char lower=cseq[k+1] >>4;
		 		sub.cseq[i]=upper | lower;
		 		sub.Ncount+= !upper + !lower;
		 		k++;
			 }
			 for(int i=0;i<len/2;i++){
					sub.Ncount+= !(cseq[k] && 0xF0) + !(cseq[k] && 0xF);
				 sub.cseq[i]	=cseq[k++];
			 }		
			 //don't forget extra half byte	 
			 if(len%2){
	 		 sub.cseq[len/2]=cseq[start/2+len/2] & 0xF0;
	 		 sub.Ncount+=!sub.cseq[len/2];
	 	 }

		 }
		}	
		return sub;
	}
	char* acgt2seq(){
		const char a[9]={'N','A','C','N','G','N','N','N','T'};
		char *seq;
		if(!(seq=(char*)malloc(seqLength+1))) exit(EXIT_FAILURE);
		for(int k=0;k<seqLength;k++){
			seq[k]=(k%2) ? a[cseq[k/2] & 0xF] : a[cseq[k/2] >>4];		 		 
		}
		return seq; 
	}
};
template <class T> class nacgt{
	//based on agct class but keeps track of positions of Ns
	//use acgt class for panel comparisons - 
	//nagct class necessary for barcode vs barcode comparisons where there may be uncertainties
	public:
	 acgt<T> A;
	 Nvector<T> N;
	 unsigned int seqLength;
	 unsigned int Ncount; 
	 nacgt(){
			seqLength=0;
			Ncount=0;
		}	
	 nacgt (char *inputSequence,int _len){
   A=acgt<T>(_len);
   N=Nvector<T>(_len);
		 seqLength=_len;
	  for (int i=0;i<seqLength;i++){
		  uint8_t k=(i%2)? 1 : 1<<4;
		  if(inputSequence[i]=='A'){
		 	 A.cseq[i/2]|=k;
		 	 continue;
		  }
		  k=k<<1;	
		  if(inputSequence[i]=='C'){
		   A.cseq[i/2]|=k;
     continue;
		  }
		  k=k<<1;	
		  if(inputSequence[i]=='G'){
		  	A.cseq[i/2]|=k;
     continue;
		  }
		  k=k<<1;
		  if(inputSequence[i]=='T'){
		  	A.cseq[i/2]|=k;
		  	continue;
		  }
    N.setBit(i);
		 }
		 Ncount= N.Ncount;
	}
	nacgt (const nacgt &a){
		Ncount=a.Ncount;
		seqLength=a.seqLength;
		A=a.A;
		N=a.N;
	}
};	
template <class T> class hamming{ 
	//functor for determining bit distances
	//class T can be uint32 or uint64 depending on the system
	//USELUT compile option exists in case there is no popcount
	
 public:
  unsigned int chunkSize=0;    //number of elements used to store the ACGT sequence
  unsigned int NchunkSize=0;   //number of elements used to store the N sequence
  unsigned int patternLength=0;//barcode size
  unsigned int ACWeight, ANWeight,NNweight;//weights for different changes
  #ifdef USELUT
   uint_fast8_t bitCountLUT[65536];
  #endif
  hamming(){
			chunkSize=0;
			NchunkSize=0;
			patternLength=0;
		 #ifdef USELUT
	  for(T i=0;i<65536;i++){
		  T temp=i;
		  bitCountLUT[i]=0;	
		 //do it the slow way
			 while(temp > 0){
			 	if (temp & 1) bitCountLUT[i]++;
			 	temp >>=  1;
			 }
			}
			#endif			
		}
  hamming(unsigned int _patternLength){
		 patternLength=_patternLength;
		 chunkSize=(patternLength%(sizeof(T)*2) ==0 )? patternLength/(sizeof(T)*2): patternLength/(sizeof(T)*2)+1;
		 NchunkSize=(patternLength%(sizeof(T)*8) ==0 )?patternLength/(sizeof(T)*8): patternLength/(sizeof(T)*8)+1;
		 //adjust chunkSize for padding if different from popcountSize
		 #ifdef USELUT
	  for(T i=0;i<65536;i++){
		  T temp=i;
		  bitCountLUT[i]=0;	
		 //do it the slow way
			 while(temp > 0){
			 	if (temp & 1) bitCountLUT[i]++;
			 	temp >>=  1;
			 }
			}
			#endif
		}
		void resize(unsigned int _patternLength){
			patternLength=_patternLength;
		 chunkSize=(patternLength%(sizeof(T)*2) ==0 )? patternLength/(sizeof(T)*2): patternLength/(sizeof(T)*2)+1;
		 NchunkSize=(patternLength%(sizeof(T)*8) ==0 )?patternLength/(sizeof(T)*8): patternLength/(sizeof(T)*8)+1;
		}	

	//the next two are the possible values for the function pointer

 unsigned int distance (nacgt<uint64_t> *a, nacgt<uint64_t> *b, unsigned int maxdist){
 	int i=0,dist=0;
 	int maxdist2=2*maxdist;
 	if(!a->A.Ncount && b->A.Ncount) dist=b->A.Ncount;     //remove these lines if the default weighting is altered
 	else if(!b->A.Ncount && a->A.Ncount) dist=a->A.Ncount;//and account for the NA transitions afterwards
			
 	//appoximate distance 
 	for(int i=0;i<chunkSize;i++){
			dist+=__builtin_popcountl(a->A.Tseq[i] ^ b->A.Tseq[i]);
			if(dist > maxdist2)return maxdist*16;
		}
		if(!a->A.Ncount || !b->A.Ncount) return dist*8;
		//rescale dist and account for the Ns
	 const int maxdist16=maxdist*16;
	 dist*=8;
	 int dNN=0;
		for(int i=0 ;i<NchunkSize;i++){
		 const int NN=__builtin_popcountl(a->N.Tseq[i] & b->N.Tseq[i]);
		 dist+=NN*15;
			if(dist > maxdist16)return maxdist16;
		 dNN+=NN;
		}
		dist+=(a->Ncount+b->Ncount-2*dNN)*4;	
		 if(dist > maxdist16)return maxdist16;
		 return dist;	
	 }
 unsigned int maxDistance (acgt<uint64_t> *panelSeqs, acgt<uint64_t> *query,int nPanelSeqs){
		//when comparing a single query with Ns we don't have to take into account the N's since
		//they give the same signal regardless of the panelSeq
		//we will initialize it anyway to get a meaningfull absolute distance
		vector<int> bestIndices(1);
		bestIndices[0]=0;
		int maxDist=query->Ncount;
 	for(int j=0;j<chunkSize;j++){
			maxDist+=__builtin_popcountl(panelSeqs[0].Tseq[j] ^ query->Tseq[j]);
		}
		for(int i=1;i<nPanelSeqs && maxDist >0 ;i++){
			int dist=query->Ncount;
		 for(int j=0;j<chunkSize && dist <= maxDist ;j++){
			 dist+=__builtin_popcountl(panelSeqs[i].Tseq[j] ^ query->Tseq[j]);
		 }
		 if(maxDist == dist){
				bestIndices.push_back(i);
			}
			else if(dist < maxDist){
				bestIndices.resize(1);
				bestIndices[0]=i;
				maxDist=dist;
			}			
		}
		//fprintf(stdout,"%.6s\n",query->acgt2seq());
		//for(int i=0;i<bestIndices.size();i++){
		//	fprintf(stdout,"%.6s %d\n",panelSeqs[bestIndices[i]].acgt2seq(),maxDist);
		//}	
		return(maxDist);	
	}
	
	int bestMatch (acgt<uint32_t> *panelSeqs, acgt<uint32_t> *query,int nPanelSeqs){
		//when comparing a single query with Ns we don't have to take into account the N's since
		//they give the same signal regardless of the panelSeq
		int bestIndex=0,nBest=1;
		int maxDist=-query->Ncount;
 	for(int j=0;j<chunkSize;j++){
			maxDist+=__builtin_popcount(panelSeqs[0].Tseq[j] ^ query->Tseq[j]);
		}
		for(int i=1;i<nPanelSeqs;i++){
			int dist=-query->Ncount;
		 for(int j=0;j<chunkSize && dist <= maxDist ;j++){
			 dist+=__builtin_popcount(panelSeqs[i].Tseq[j] ^ query->Tseq[j]);
		 }
		 if(maxDist == dist){
				if(!dist)return -1; //dupe found
				nBest++;
			}
			else if(dist < maxDist){
				nBest=1;
				bestIndex=i;
				maxDist=dist;
			}			
		}
		
		if(maxDist ==0 && nBest==1){
			return bestIndex;
		}
		return -1;	
	}	
	unsigned int bestMatch (vector<string> &panelSeqs, char *query,int nPanelSeqs){
		//when comparing a single query with Ns we don't have to take into account the N's since
		//they give the same signal regardless of the panelSeq
		//we will initialize it anyway to get a meaningfull absolute distance
		int bestIndex=0,nBest=1;
		int maxDist=0;
  for (unsigned int i = 0; i < patternLength; i++){
	 	if(query[i] != 'N' && panelSeqs[0][i] !=query[i]){
    maxDist++;
	 	}	
	 }
		for(int i=1;i<nPanelSeqs ;i++){
   int dist=0;
   for (unsigned int j = 0; j < patternLength && dist<=maxDist; j++){
				if(query[j] != 'N' && panelSeqs[i][j] !=query[j]){
     dist++;
				}	
			}	
		 if(maxDist == dist){				
				//if(!dist)return -1; //dupe found
				nBest++;
			}
			else if(dist < maxDist){
				nBest=1;
				bestIndex=i;
				maxDist=dist;
			}			
		}

		if(maxDist ==0 && nBest==1){

			return bestIndex;
		}
		return -1;	
	}
		unsigned int bestMatch (vector<string> &panelSeqs, char *query,int nPanelSeqs,int mismatchTol,int NTol){
		//when comparing a single query with Ns we don't have to take into account the N's since
		//they give the same signal regardless of the panelSeq
		//we will initialize it anyway to get a meaningfull absolute distance
		int bestIndex=0,nBest=1;
		int maxDist=0;
		int nNs=0;
  for (unsigned int i = 0; i < patternLength; i++){
	 	if(query[i] != 'N' && panelSeqs[0][i] !=query[i]){
    maxDist++;
	 	}
	 	else if (query[i] == 'N') nNs++;
	 	if(nNs > NTol)return -1;
	 }
		for(int i=1;i<nPanelSeqs ;i++){
   int dist=0;
   for (unsigned int j = 0; j < patternLength && dist<=maxDist; j++){
				if(query[j] != 'N' && panelSeqs[i][j] !=query[j]){
     dist++;
				}	
			}	
		 if(maxDist == dist){				
				//if(!dist)return -1; //dupe found
				nBest++;
			}
			else if(dist < maxDist){
				nBest=1;
				bestIndex=i;
				maxDist=dist;
			}			
		}
		if(maxDist <=mismatchTol && nBest==1){
			return bestIndex;
		}
		return -1;	
	}
	 unsigned int distance (nacgt<uint32_t> *a, nacgt<uint32_t> *b, unsigned int maxdist){
 	int i=0,dist=0;
 	int maxdist2=2*maxdist;
 	//appoximate distance 
 	for(int i=0;i<chunkSize;i++){
			dist+=__builtin_popcount(a->A.Tseq[i] ^ b->A.Tseq[i]);
			if(dist > maxdist2)return maxdist*16;
		}
		//rescale dist and account for the Ns
	 const int maxdist16=maxdist*16;
	 dist*=8;
	 int dNN=0;
		for(int i=0 ;i<NchunkSize;i++){
			const int NN=__builtin_popcount(a->N.Tseq[i] & b->N.Tseq[i]);
			dist+=NN*15;
			if(dist > maxdist16)return maxdist16;
			dNN+=NN;
		}
		dist+=(a->Ncount+b->Ncount-2*dNN)*4;	
		if(dist > maxdist16)return maxdist16;
		return dist;	
	}

 unsigned int slowDistance(const char *a, const char *b, unsigned int len,int maxdist){
	 int dist=0;
	 const int maxdist16=16*maxdist;			
  for (unsigned int i = 0; i < len && dist<maxdist16; i++){

	 	if(a[i] !=b[i]){
	 		if(a[i] == 'N' || b[i] == 'N'){
					dist+=12;
				}
	 		else{
				 dist+=16;
				}
	 	}
	 	else if(a[i] == 'N'){
				dist+=15;
			}		
	 }
	 if(dist > maxdist16) return maxdist16;
  return(dist);
	}

	#ifdef USELUT
	unsigned int hpopcount(nacgt<uint32_t> *a, nacgt<uint32_t> *b, unsigned int maxdist){
 	int i=0,dist=0; 
 	int maxdist2=2*maxdist;
 	if(!a->A.Ncount && b->A.Ncount) dist=b->A.Ncount;
 	else if(!b->A.Ncount && a->A.Ncount) dist=a->A.Ncount;
			
 	//appoximate distance 
 	for(int i=0;i<chunkSize;i++){
			dist+=__builtin_popcountl(a->A.Tseq[i] ^ b->A.Tseq[i]);
			if(dist > maxdist2)return maxdist*16;
		}
		if(!a->A.Ncount || !b->A.Ncount) return dist*8;
		
		//rescale dist and account for the Ns
	 const int maxdist16=maxdist*16;
	 dist*=8;
	 int dNN=0;
		for(int i=0 ;i<NchunkSize;i++){
			uint32_t diff= a->N.Tseq[i] & b->N.Tseq[i];
			const int NN=bitCountLUT[diff&0xFFFF] + bitCountLUT[diff >> 16] ;
			dist+=NN*15;
			if(dist > maxdist16)return maxdist16;
			dNN+=NN;
		}
		dist+=(a->Ncount+b->Ncount-2*dNN)*4;	
		if(dist > maxdist16)return maxdist16;
		return dist;	
	}
	#endif	
};	

static void printBits(void *ptr, char nBytes){
	
	char *cPtr = (char*) ptr;
	for(int i=0;i<nBytes;i++){
		//for example to gt second byte
		//1100 0000 >> 6 = 0000 0011
		//               & 0000 0001     
		//                 0000 0001
		for(int b=7;b>=0;b--){
			printf("%d",(cPtr[i] >> b) & 1);
		}	
		printf(" ");
	}
	printf("\n");
}
static char* reverseBytes(void *ptr, char nBytes){
	char *cptr=(char*) malloc(nBytes);
	memcpy(cptr,ptr,nBytes);
	char *left=(char*) cptr;
	char *right=left+nBytes-1;
	while(left <right){
		char temp=*left;
		*left++=*right;
		*right--=temp;
	}
	return(cptr);
}
void printBitsEndian(void *ptr, char nBytes, char endian){
	if(nBytes == 1 || !endian) {
		printBits(ptr,nBytes);
	 return;
	}
	
	uint16_t test=1;
	//0000 0000  0000 0001 BigEndian
	//0000 0001  0000 0000 Little Endian
	char *a = (char*) & test; //pointer cast to byte type
	char BigEndianSystem=a[1] & 1; //look at second byte
	if(endian == -1 && BigEndianSystem) printBits(ptr,nBytes);
	else if(endian ==-1 && !BigEndianSystem){
		char *revptr=reverseBytes(ptr,nBytes);
		printBits(revptr,nBytes);
		free(revptr);
	}
	else if(endian=1 && BigEndianSystem){
		char *revptr=reverseBytes(ptr,nBytes);
		reverseBytes(ptr,nBytes);
		printBits(revptr,nBytes);
		free(revptr);
	}
	else{
		printBits(ptr,nBytes);
	}	
}


