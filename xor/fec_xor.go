package main

/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <immintrin.h>

// 1		: XOR operation with SIMD
// Other 	: XOR operation without SIMD

#define SIMD 1

char** xorENC(char** src, int k, int s);
void vector_xor(char *dst, char *src1, char *src2, int iter_num);
void xorDEC(char** enc, int idx, int k, int s);
int sim(char** enc, int k, int s);
void print(char** arr, int k, int s);
char** make2Darray(int row, int col);
void fileWrite(char** enc, int cnt, int col, FILE *wp);
int fileRead(char** data, int row, int col, FILE *rp);
void freeArray(char** enc, int row);

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

char** make2Darray(int row, int col)
{
	char** data = (char**)malloc(row * sizeof(char*)) ;

	for(int i = 0 ; i < row ; i++)
		data[i] = (char*)malloc(col * sizeof(char)) ;

	return data ;
}

void fileWrite(char** enc, int cnt, int col, FILE *wp)
{
	int max_r = (int)ceil((float)cnt / col) ;

	for(int i = 0 ; i < max_r ; i++)
		fwrite(enc[i], col, 1, wp) ;
}

int fileRead(char** data, int row, int col, FILE *rp)
{
	int cnt = fread(*data, 1, row * col, rp) ;
	return cnt ;
}

void freeArray(char** enc, int row)
{
	for (int i = 0 ; i < row + 1 ; i++)
				free(enc[i]) ;

	free(enc) ;
}

*/
import "C"

import (
	"flag"
	"fmt"
	"os"
)

// 32KIB
var row = flag.Int("r", 32, "Size of row")
var col = flag.Int("c", 1024, "Size of column")

func init() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage : %s [-flags] inputFile\n\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "Valid flags:\n")
		flag.PrintDefaults()
	}
}

func main() {
	flag.Parse()
	args := flag.Args()
	if len(args) != 1 {
		fmt.Fprintf(os.Stderr, "Error: No input filename given\n")
		flag.Usage()
		os.Exit(1)
	}
	fname := args[0]

	var enc, data **C.char
	var idx, cnt C.int
	var rp, wp *C.FILE

	data = C.make2Darray(C.int(*row), C.int(*col))

	rp = C.fopen(C.CString(fname), C.CString("rb"))
	wp = C.fopen(C.CString("result.txt"), C.CString("wb"))

	for C.feof(rp) == 0 {
		cnt = C.fileRead(data, C.int(*row), C.int(*col), rp)

		enc = C.xorENC(data, C.int(*row), C.int(*col))

		idx = C.sim(enc, C.int(*row), C.int(*col))

		C.xorDEC(enc, idx, C.int(*row), C.int(*col))

		C.fileWrite(enc, cnt, C.int(*col), wp)

		C.freeArray(enc, C.int(*row))
	}
}
