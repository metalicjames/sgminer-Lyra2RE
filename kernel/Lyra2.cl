

/*Blake2b IV Array*/
__constant sph_u64 blake2b_IV[8] =
{
  0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL,
  0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
  0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
  0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL
};

/*Blake2b's rotation*/
static inline sph_u64 rotr64( const sph_u64 w, const unsigned c ){
    return ( w >> c ) | ( w << ( 64 - c ) );
}

/*Blake2b's G function*/
#define G(a,b,c,d) \
  do { \
a += b; d ^= a; d = SPH_ROTR64(d, 32); \
c += d; b ^= c; b = SPH_ROTR64(b, 24); \
a += b; d ^= a; d = SPH_ROTR64(d, 16); \
c += d; b ^= c; b = SPH_ROTR64(b, 63); \
  } while(0)


/*One Round of the Blake2b's compression function*/
#define round_lyra(v)  \
 do { \
    G(v[ 0],v[ 4],v[ 8],v[12]); \
    G(v[ 1],v[ 5],v[ 9],v[13]); \
    G(v[ 2],v[ 6],v[10],v[14]); \
    G(v[ 3],v[ 7],v[11],v[15]); \
    G(v[ 0],v[ 5],v[10],v[15]); \
    G(v[ 1],v[ 6],v[11],v[12]); \
    G(v[ 2],v[ 7],v[ 8],v[13]); \
    G(v[ 3],v[ 4],v[ 9],v[14]); \
 } while(0)


#define reduceDuplexRowSetup(rowIn, rowInOut, rowOut) \
 do { \
	for (int i = 0; i < 8; i++) \
	{ \
\
		for (int j = 0; j < 12; j++) {state[j] ^= Matrix[rowIn][12 * i + j] + Matrix[rowInOut][12 * i + j];} \
		round_lyra(state); \
		for (int j = 0; j < 12; j++) {Matrix[rowOut][j + 84 - 12 * i] = Matrix[rowIn][12 * i + j] ^ state[j];} \
\
		Matrix[rowInOut][0 + 12 * i] ^= state[11]; \
		Matrix[rowInOut][1 + 12 * i] ^= state[0]; \
		Matrix[rowInOut][2 + 12 * i] ^= state[1]; \
		Matrix[rowInOut][3 + 12 * i] ^= state[2]; \
		Matrix[rowInOut][4 + 12 * i] ^= state[3]; \
		Matrix[rowInOut][5 + 12 * i] ^= state[4]; \
		Matrix[rowInOut][6 + 12 * i] ^= state[5]; \
		Matrix[rowInOut][7 + 12 * i] ^= state[6]; \
		Matrix[rowInOut][8 + 12 * i] ^= state[7]; \
		Matrix[rowInOut][9 + 12 * i] ^= state[8]; \
		Matrix[rowInOut][10 + 12 * i] ^= state[9]; \
		Matrix[rowInOut][11 + 12 * i] ^= state[10]; \
	} \
 \
 } while (0)

 #define reduceDuplexRow(rowIn, rowInOut, rowOut) \
 do { \
 \
	 for (int i = 0; i < 8; i++) \
	 	 { \
		 for (int j = 0; j < 12; j++) \
			 state[j] ^= Matrix[rowIn][12 * i + j] + Matrix[rowInOut][12 * i + j]; \
 \
		 round_lyra(state); \
		 for (int j = 0; j < 12; j++) Matrix[rowOut][j + 12 * i] ^= state[j]; \
		 Matrix[rowInOut][0 + 12 * i] ^= state[11]; \
		 Matrix[rowInOut][1 + 12 * i] ^= state[0]; \
		 Matrix[rowInOut][2 + 12 * i] ^= state[1]; \
		 Matrix[rowInOut][3 + 12 * i] ^= state[2]; \
		 Matrix[rowInOut][4 + 12 * i] ^= state[3]; \
		 Matrix[rowInOut][5 + 12 * i] ^= state[4]; \
		 Matrix[rowInOut][6 + 12 * i] ^= state[5]; \
		 Matrix[rowInOut][7 + 12 * i] ^= state[6]; \
		 Matrix[rowInOut][8 + 12 * i] ^= state[7]; \
		 Matrix[rowInOut][9 + 12 * i] ^= state[8]; \
		 Matrix[rowInOut][10 + 12 * i] ^= state[9]; \
		 Matrix[rowInOut][11 + 12 * i] ^= state[10]; \
	 	 } \
 \
 } while (0)

#define absorbblock(in) do { \
	state[0] ^= Matrix[in][0]; \
	state[1] ^= Matrix[in][1]; \
	state[2] ^= Matrix[in][2]; \
	state[3] ^= Matrix[in][3]; \
	state[4] ^= Matrix[in][4]; \
	state[5] ^= Matrix[in][5]; \
	state[6] ^= Matrix[in][6]; \
	state[7] ^= Matrix[in][7]; \
	state[8] ^= Matrix[in][8]; \
	state[9] ^= Matrix[in][9]; \
	state[10] ^= Matrix[in][10]; \
	state[11] ^= Matrix[in][11]; \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
} while(0)
