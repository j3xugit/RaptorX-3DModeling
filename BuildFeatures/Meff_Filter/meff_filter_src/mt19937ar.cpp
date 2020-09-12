#include "mt19937ar.h"

// Period parameters //
//#define N 624
//#define M 397
//#define MATRIX_A 0x9908b0dfUL   // constant vector a //
//#define UPPER_MASK 0x80000000UL // most significant w-r bits //
//#define LOWER_MASK 0x7fffffffUL // least significant r bits //
//static unsigned long mt[N]; // the array for the state vector  //
//static int mti=N+1; // mti==N+1 means mt[N] is not initialized //

mt19937::mt19937(void)
{
	N=624;
	M=397;
	MATRIX_A=0x9908b0dfUL;
	UPPER_MASK=0x80000000UL;
	LOWER_MASK=0x7fffffffUL;
	mt=new unsigned long[N];
	mti=N+1;

	idum=0;   //__091210__//
	ival=0.0; //__091210__//
}
mt19937::~mt19937(void)
{
}

//--------------- initialize ---------------//
// initializes mt[N] with a seed //
void mt19937::init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. //
        // In the previous versions, MSBs of the seed affect   //
        // only MSBs of the array mt[].                        //
        // 2002/01/09 modified by Makoto Matsumoto             //
        mt[mti] &= 0xffffffffUL;
        // for >32 bit machines //
    }
}


//-----------------------------------------------------------------
//    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
//    init_by_array(init, length);
//------------------------------------------------------------------//
// initialize by an array with array-length //
// init_key is the array for initializing keys //
// key_length is its length //
// slight change for C++, 2004/2/26 //
void mt19937::init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; // non linear //
        mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines //
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; // non linear //
        mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines //
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }
    mt[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array //
}

//------------------ random number generate -------------------//
//-- generates a random number on [0,0xffffffff]-interval --//
unsigned long mt19937::genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1 //

    if (mti >= N) { // generate N words at one time //
        int kk;

        if (mti == N+1)   // if init_genrand() has not been called, //
            init_genrand(5489UL); // a default initial seed is used //

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    // Tempering //
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

//-- generates a random number on [0,0x7fffffff]-interval --//
long mt19937::genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

//-- generates a random number on [0,1]-real-interval --//
double mt19937::genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    // divided by 2^32-1 // 
}

//-- generates a random number on [0,1)-real-interval --//
double mt19937::genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    // divided by 2^32 //
}

//-- generates a random number on (0,1)-real-interval --//
double mt19937::genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    // divided by 2^32 //
}

//-- generates a random number on [0,1) with 53-bit resolution --//
double mt19937::genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
}


//---------------- additional -------------------//
//[input] must positive
int mt19937::choose_rand(double *input,int total)
{
	int i;
	double value;
	double totval;
	totval=0.0;
	for(i=0;i<total;i++)totval+=input[i];
	value=genrand_real3()*totval;
	totval=0.0;
	for(i=0;i<total;i++)
	{
		totval+=input[i];
		if(totval>value)return i;
	}
	return -1;
}
//sampling from N(0,1) //Box-Muller method
double mt19937::normal_distribution(void)
{
	double v1,v2;
	double rsq,fac;

	if(idum==0)
	{
		for(;;)
		{
			v1=2.0*genrand_real3()-1.0;
			v2=2.0*genrand_real3()-1.0;
			rsq=v1*v1+v2*v2;
			if(rsq>=1.0||rsq==0.0)continue;
			else break;
		}
		fac=sqrt(-2.0*log(rsq)/rsq);
		ival=v1*fac;
		idum=1;
		return v2*fac;
	}
	else
	{
		idum=0;
		return ival;
	}
}
