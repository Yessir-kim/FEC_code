#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <immintrin.h>

// 1		: XOR operation with SIMD
// Other 	: XOR operation without SIMD

#define SIMD 1 

#define LOG 0 

char** xorENC(char** src, int k, int s);
void vector_xor(char *dst, char *src1, char *src2, int iter_num);
void xorDEC(char** enc, int idx, int k, int s);
int sim(char** enc, int k, int s);
void print(char** arr, int k, int s);

int main(int argc, char *argv[])
{
	if(argc < 4) 
	{
		printf("Usage: %s <row> <col> <fileName>\n", argv[0]) ;
		return 0 ;
	}

	// row * col = block size 
	int row = atoi(argv[1]) ;
	int col = atoi(argv[2]) ;
	char* fileName = argv[3] ;
	char** enc ; 
	int idx, cnt, max_r ;
	FILE *rp, *wp ; 

	char** data = (char**)malloc(row * sizeof(char*)) ;

	for(int i = 0 ; i < row ; i++)
		data[i] = (char*)malloc(col * sizeof(char)) ;

	// GET FILE POINTER (read & write)
	rp = fopen(fileName, "rb") ;
	wp = fopen("result.txt", "wb") ;

	while(!feof(rp))
	{
		// LOAD DATA FROM FILE
		cnt = fread(*data, 1, row * col, rp) ; 

		#if LOG == 1

		printf("--------- LOAD ---------\n");
		print(data, row, col) ;

		#endif

		// ENCODING
		enc = xorENC(data, row, col) ;
	
		#if LOG == 1

		printf("--------- ENCODE ---------\n");
		print(enc, row + 1, col) ;
	
		#endif

		// LOSS 
		idx = sim(enc, row, col) ;
	
		#if LOG == 1

		printf("--------- LOSS(%d) ---------\n", idx);
		print(enc, row, col) ;
	
		#endif
	
		// DECODING
		xorDEC(enc, idx, row, col) ;
	
		#if LOG == 1

		printf("--------- DECODE ---------\n");
		print(enc, row, col) ;

		#endif

		max_r = (int)ceil((float)cnt / col) ;

		// WRITE DATA
		for(int i = 0 ; i < max_r ; i++)
			fwrite(enc[i], col, 1, wp) ;
		
		for (int i = 0 ; i < row + 1 ; i++)
        	free(enc[i]) ;

		free(enc) ;
	}

	fclose(wp) ;
	fclose(rp) ;

	for (int i = 0 ; i < row ; i++) 
		free(data[i]) ;
	
	free(data) ;
}

char** xorENC(char** src, int k, int s)
{
	clock_t start, end ;
	
	char** res = (char**)malloc((k + 1) * sizeof(char*));

	for(int i = 0 ; i < k + 1 ; i++)
	{
		res[i] = (char*)malloc(s * sizeof(char)) ;
		memset(res[i], '\0', s * sizeof(char)) ;
	}

	start = clock() ;

	#if SIMD == 1
	
	memcpy(res[k], src[0], s) ; // init

	for(int i = 0 ; i < k ; i++)
	{
		memcpy(res[i], src[i], s) ; // to copy original data to result array
		
		if(i != 0) vector_xor(res[k], res[k], src[i], s) ; // XOR operation
	}

	#else
	
	for(int i = 0 ; i < k ; i++)
	{
		for(int j = 0 ; j < s ; j++)
		{
			res[i][j] = src[i][j] ; // to copy original data to result array

			if(i == 0) res[k][j] = src[i][j] ; // init
			else res[k][j] = res[k][j] ^ src[i][j] ; // XOR operation
		}
	
	}

	#endif

	end = clock() ;

	printf("xorENC execution time, %f\n", (double)(end - start)/CLOCKS_PER_SEC) ;

	return res ;
}

void xorDEC(char** enc, int idx, int k, int s)
{
	clock_t start, end ;

	start = clock() ;

	#if SIMD == 1
	
	memcpy(enc[idx], enc[k], s) ; // init 

	for(int i = 0 ; i < k ; i++)
	{
		if(i != idx) 
			vector_xor(enc[idx], enc[idx], enc[i], s) ; // XOR operation
	}

	#else
	
	for(int i = 0 ; i < k ; i++)
	{
		for(int j = 0 ; j < s ; j++)
			if(i != idx) enc[k][j] = enc[k][j] ^ enc[i][j] ; // XOR operation
	}

	for(int i = 0 ; i < s ; i++)
		enc[idx][i] = enc[k][i] ; // restore data

	#endif

	end = clock() ;

	printf("xorDEC execution time, %f\n", (double)(end - start)/CLOCKS_PER_SEC) ;
}

int sim(char** enc, int k, int s)
{
	srand(time(NULL)) ;
	int rnd = rand() % k ; // loss idx 

	for(int i = 0 ; i < s ; i++)
		enc[rnd][i] = '0' ; // fill with the zero value

	return rnd ;
}

void print(char** arr, int k, int s)
{
	for(int i = 0 ; i < k ; i++)
    { 
		for(int j = 0 ; j < s ; j++)
			printf("%c ", arr[i][j]) ;
		
		printf("\n") ;
	}
}

void vector_xor(char *dst, char *src1, char *src2, int iter_num)
{
	__m128i a, b ;
	__m128i result ;

	for (int i = 0 ; i < iter_num ; i += 16) 
	{
		a = _mm_loadu_si128((__m128i const*)(src1+i)) ;
		b = _mm_loadu_si128((__m128i const*)(src2+i)) ;
		result = _mm_xor_si128(a, b) ; 
		_mm_storeu_si128((__m128i_u *)(dst+i), result) ;
	}
}
