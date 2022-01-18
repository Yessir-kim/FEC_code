package main

/*
#define GF_BITS  8

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include "rs.h"

#ifdef  TEST
#define DEB(x)
#define DDB(x) x
#define DEBUG   0

#include <sys/time.h>
#define DIFF_T(a,b) \
    (1+ 1000000*(a.tv_sec - b.tv_sec) + (a.tv_usec - b.tv_usec) )

#define TICK(t) \
    {struct timeval x ; \
    gettimeofday(&x, NULL) ; \
    t = x.tv_usec + 1000000* (x.tv_sec & 0xff ) ; \
    }
#define TOCK(t) \
    { u_long t1 ; TICK(t1) ; \
      if (t1 < t) t = 256000000 + t1 - t ; \
      else t = t1 - t ; \
      if (t == 0) t = 1 ;}

u_long ticks[10];
#else
#define DEB(x)
#define DDB(x)
#define TICK(x)
#define TOCK(x)
#endif

#if (GF_BITS != 8)
#error "GF_BITS must be 8"
#endif
typedef unsigned char gf;

#define GF_SIZE ((1 << GF_BITS) - 1)

static char *allPp[] = {
    NULL,
    NULL,
    "111",
    "1101",
    "11001",
    "101001",
    "1100001",
    "10010001",
    "101110001",
    "1000100001",
    "10010000001",
    "101000000001",
    "1100101000001",
    "11011000000001",
    "110000100010001",
    "1100000000000001",
    "11010000000010001"
};


static gf gf_exp[2*GF_SIZE];
static int gf_log[GF_SIZE + 1];
static gf inverse[GF_SIZE+1];

static inline gf
modnn(int x)
{
    while (x >= GF_SIZE) {
    x -= GF_SIZE;
    x = (x >> GF_BITS) + (x & GF_SIZE);
    }
    return x;
}

#define SWAP(a,b,t) {t tmp; tmp=a; a=b; b=tmp;}

static gf gf_mul_table[(GF_SIZE + 1)*(GF_SIZE + 1)]
#ifdef WINDOWS
__attribute__((aligned (16)))
#else
__attribute__((aligned (256)))
#endif
;

#define gf_mul(x,y) gf_mul_table[(x<<8)+y]

#define USE_GF_MULC register gf * __gf_mulc_
#define GF_MULC0(c) __gf_mulc_ = &gf_mul_table[(c)<<8]
#define GF_ADDMULC(dst, x) dst ^= __gf_mulc_[x]
#define GF_MULC(dst, x) dst = __gf_mulc_[x]

static void
init_mul_table(void)
{
    int i, j;
    for (i=0; i< GF_SIZE+1; i++)
    for (j=0; j< GF_SIZE+1; j++)
        gf_mul_table[(i<<8)+j] = gf_exp[modnn(gf_log[i] + gf_log[j]) ] ;

    for (j=0; j< GF_SIZE+1; j++)
    gf_mul_table[j] = gf_mul_table[j<<8] = 0;
}

static void
generate_gf(void)
{
    int i;
    gf mask;
    char *Pp =  allPp[GF_BITS] ;

    mask = 1;
    gf_exp[GF_BITS] = 0;

    for (i = 0; i < GF_BITS; i++, mask <<= 1 ) {
    gf_exp[i] = mask;
    gf_log[gf_exp[i]] = i;

    if ( Pp[i] == '1' )
        gf_exp[GF_BITS] ^= mask;
    }

    gf_log[gf_exp[GF_BITS]] = GF_BITS;

    mask = 1 << (GF_BITS - 1 ) ;
    for (i = GF_BITS + 1; i < GF_SIZE; i++) {
    if (gf_exp[i - 1] >= mask)
        gf_exp[i] = gf_exp[GF_BITS] ^ ((gf_exp[i - 1] ^ mask) << 1);
    else
        gf_exp[i] = gf_exp[i - 1] << 1;
    gf_log[gf_exp[i]] = i;
    }

    gf_log[0] = GF_SIZE ;

    for (i = 0 ; i < GF_SIZE ; i++)
    gf_exp[i + GF_SIZE] = gf_exp[i] ;

    inverse[0] = 0 ;
    inverse[1] = 1;
    for (i=2; i<=GF_SIZE; i++)
    inverse[i] = gf_exp[GF_SIZE-gf_log[i]];
}

#if 0
#define addmul(dst, src, c, sz) \
    if (c != 0) addmul1(dst, src, c, sz)
#endif



#define UNROLL 16
static void
slow_addmul1(gf *dst1, gf *src1, gf c, int sz)
{
    USE_GF_MULC ;
    register gf *dst = dst1, *src = src1 ;
    gf *lim = &dst[sz - UNROLL + 1] ;

    GF_MULC0(c) ;

#if (UNROLL > 1)
    for (; dst < lim ; dst += UNROLL, src += UNROLL ) {
    GF_ADDMULC( dst[0] , src[0] );
    GF_ADDMULC( dst[1] , src[1] );
    GF_ADDMULC( dst[2] , src[2] );
    GF_ADDMULC( dst[3] , src[3] );
#if (UNROLL > 4)
    GF_ADDMULC( dst[4] , src[4] );
    GF_ADDMULC( dst[5] , src[5] );
    GF_ADDMULC( dst[6] , src[6] );
    GF_ADDMULC( dst[7] , src[7] );
#endif
#if (UNROLL > 8)
    GF_ADDMULC( dst[8] , src[8] );
    GF_ADDMULC( dst[9] , src[9] );
    GF_ADDMULC( dst[10] , src[10] );
    GF_ADDMULC( dst[11] , src[11] );
    GF_ADDMULC( dst[12] , src[12] );
    GF_ADDMULC( dst[13] , src[13] );
    GF_ADDMULC( dst[14] , src[14] );
    GF_ADDMULC( dst[15] , src[15] );
#endif
    }
#endif
    lim += UNROLL - 1 ;
    for (; dst < lim; dst++, src++ )
    GF_ADDMULC( *dst , *src );
}

# define addmul1 slow_addmul1

static void addmul(gf *dst, gf *src, gf c, int sz) {
    if (c != 0) addmul1(dst, src, c, sz);
}


#if 0
#define mul(dst, src, c, sz) \
    do { if (c != 0) mul1(dst, src, c, sz); else memset(dst, 0, c); } while(0)
#endif

#define UNROLL 16
static void
slow_mul1(gf *dst1, gf *src1, gf c, int sz)
{
    USE_GF_MULC ;
    register gf *dst = dst1, *src = src1 ;
    gf *lim = &dst[sz - UNROLL + 1] ;

    GF_MULC0(c) ;

#if (UNROLL > 1)
    for (; dst < lim ; dst += UNROLL, src += UNROLL ) {
    GF_MULC( dst[0] , src[0] );
    GF_MULC( dst[1] , src[1] );
    GF_MULC( dst[2] , src[2] );
    GF_MULC( dst[3] , src[3] );
#if (UNROLL > 4)
    GF_MULC( dst[4] , src[4] );
    GF_MULC( dst[5] , src[5] );
    GF_MULC( dst[6] , src[6] );
    GF_MULC( dst[7] , src[7] );
#endif
#if (UNROLL > 8)
    GF_MULC( dst[8] , src[8] );
    GF_MULC( dst[9] , src[9] );
    GF_MULC( dst[10] , src[10] );
    GF_MULC( dst[11] , src[11] );
    GF_MULC( dst[12] , src[12] );
    GF_MULC( dst[13] , src[13] );
    GF_MULC( dst[14] , src[14] );
    GF_MULC( dst[15] , src[15] );
#endif
    }
#endif
    lim += UNROLL - 1 ;
    for (; dst < lim; dst++, src++ )
    GF_MULC( *dst , *src );
}

# define mul1 slow_mul1

static inline void mul(gf *dst, gf *src, gf c, int sz) {
    if (c != 0) mul1(dst, src, c, sz); else memset(dst, 0, c);
}

DEB( int pivloops=0; int pivswaps=0 ;)
    static int
invert_mat(gf *src, int k)
{
    gf c, *p ;
    int irow, icol, row, col, i, ix ;

    int error = 1 ;
    int indxc[k];
    int indxr[k];
    int ipiv[k];
    gf id_row[k];

    memset(id_row, 0, k*sizeof(gf));
    DEB( pivloops=0; pivswaps=0 ; )

    for (i = 0; i < k ; i++)
        ipiv[i] = 0 ;

    for (col = 0; col < k ; col++) {
    gf *pivot_row ;

    irow = icol = -1 ;
    if (ipiv[col] != 1 && src[col*k + col] != 0) {
        irow = col ;
        icol = col ;
        goto found_piv ;
    }
    for (row = 0 ; row < k ; row++) {
        if (ipiv[row] != 1) {
        for (ix = 0 ; ix < k ; ix++) {
            DEB( pivloops++ ; )
            if (ipiv[ix] == 0) {
                if (src[row*k + ix] != 0) {
                irow = row ;
                icol = ix ;
                goto found_piv ;
                }
            } else if (ipiv[ix] > 1) {
                fprintf(stderr, "singular matrix\n");
                goto fail ;
            }
        }
        }
    }
    if (icol == -1) {
        fprintf(stderr, "XXX pivot not found!\n");
        goto fail ;
    }
 found_piv:
    ++(ipiv[icol]) ;

    if (irow != icol) {
        for (ix = 0 ; ix < k ; ix++ ) {
        SWAP( src[irow*k + ix], src[icol*k + ix], gf) ;
        }
    }
    indxr[col] = irow ;
    indxc[col] = icol ;
    pivot_row = &src[icol*k] ;
    c = pivot_row[icol] ;
    if (c == 0) {
        fprintf(stderr, "singular matrix 2\n");
        goto fail ;
    }
    if (c != 1 ) {

        DEB( pivswaps++ ; )
        c = inverse[ c ] ;
        pivot_row[icol] = 1 ;
        for (ix = 0 ; ix < k ; ix++ )
        pivot_row[ix] = gf_mul(c, pivot_row[ix] );
    }

    id_row[icol] = 1;
    if (memcmp(pivot_row, id_row, k*sizeof(gf)) != 0) {
        for (p = src, ix = 0 ; ix < k ; ix++, p += k ) {
        if (ix != icol) {
            c = p[icol] ;
            p[icol] = 0 ;
            addmul(p, pivot_row, c, k );
        }
        }
    }
    id_row[icol] = 0;
    }
    for (col = k-1 ; col >= 0 ; col-- ) {
    if (indxr[col] <0 || indxr[col] >= k)
        fprintf(stderr, "AARGH, indxr[col] %d\n", indxr[col]);
    else if (indxc[col] <0 || indxc[col] >= k)
        fprintf(stderr, "AARGH, indxc[col] %d\n", indxc[col]);
    else
        if (indxr[col] != indxc[col] ) {
        for (row = 0 ; row < k ; row++ ) {
            SWAP( src[row*k + indxr[col]], src[row*k + indxc[col]], gf) ;
        }
        }
    }
    error = 0 ;
 fail:
    return error ;
}

static int fec_initialized = 0 ;

void fec_init(void)
{
    TICK(ticks[0]);
    generate_gf();
    TOCK(ticks[0]);
    DDB(fprintf(stderr, "generate_gf took %ldus\n", ticks[0]);)
    TICK(ticks[0]);
    init_mul_table();
    TOCK(ticks[0]);
    DDB(fprintf(stderr, "init_mul_table took %ldus\n", ticks[0]);)
    fec_initialized = 1 ;
}


#ifdef PROFILE
static long long rdtsc(void)
{
    unsigned long low, hi;
    asm volatile ("rdtsc" : "=d" (hi), "=a" (low));
    return ( (((long long)hi) << 32) | ((long long) low));
}

void print_matrix1(gf* matrix, int nrows, int ncols) {
    int i, j;
    printf("matrix (%d,%d):\n", nrows, ncols);
    for(i = 0; i < nrows; i++) {
        for(j = 0; j < ncols; j++) {
            printf("%6d ", matrix[i*ncols + j]);
        }
        printf("\n");
    }
}

void print_matrix2(gf** matrix, int nrows, int ncols) {
    int i, j;
    printf("matrix (%d,%d):\n", nrows, ncols);
    for(i = 0; i < nrows; i++) {
        for(j = 0; j < ncols; j++) {
            printf("%6d ", matrix[i][j]);
        }
        printf("\n");
    }
}

#endif


static gf galExp(gf a, gf n) {
    int logA;
    int logResult;
    if(0 == n) {
        return 1;
    }
    if(0 == a) {
        return 0;
    }
    logA = gf_log[a];
    logResult = logA * n;
    while(logResult >= 255) {
        logResult -= 255;
    }

    return gf_exp[logResult];
}

static inline gf galMultiply(gf a, gf b) {
    return gf_mul_table[ ((int)a << 8) + (int)b ];
}

static gf* vandermonde(int nrows, int ncols) {
    int row, col, ptr;
    gf* matrix = (gf*)RS_MALLOC(nrows * ncols);
    if(NULL != matrix) {
        ptr = 0;
        for(row = 0; row < nrows; row++) {
            for(col = 0; col < ncols; col++) {
                matrix[ptr++] = galExp((gf)row, (gf)col);
            }
        }
    }

    return matrix;
}


static gf* sub_matrix(gf* matrix, int rmin, int cmin, int rmax, int cmax,  int nrows, int ncols) {
    int i, j, ptr = 0;
    gf* new_m = (gf*)RS_MALLOC( (rmax-rmin) * (cmax-cmin) );
    if(NULL != new_m) {
        for(i = rmin; i < rmax; i++) {
            for(j = cmin; j < cmax; j++) {
                new_m[ptr++] = matrix[i*ncols + j];
            }
        }
    }

    return new_m;
}


static gf* multiply1(gf *a, int ar, int ac, gf *b, int br, int bc) {
    gf *new_m, tg;
    int r, c, i, ptr = 0;

    assert(ac == br);
    new_m = (gf*)RS_CALLOC(1, ar*bc);
    if(NULL != new_m) {

        for(r = 0; r < ar; r++) {
            for(c = 0; c < bc; c++) {
                tg = 0;
                for(i = 0; i < ac; i++) {
                    tg ^= galMultiply(a[r*ac+i], b[i*bc+c]);
                }

                new_m[ptr++] = tg;
            }
        }

    }

    return new_m;
}

static inline int code_some_shards(gf* matrixRows, gf** inputs, gf** outputs,
        int dataShards, int outputCount, int byteCount) {
    gf* in;
    int iRow, c;
    for(c = 0; c < dataShards; c++) {
        in = inputs[c];
        for(iRow = 0; iRow < outputCount; iRow++) {
            if(0 == c) {
                mul(outputs[iRow], in, matrixRows[iRow*dataShards+c], byteCount);
            } else {
                addmul(outputs[iRow], in, matrixRows[iRow*dataShards+c], byteCount);
            }
        }
    }

    return 0;
}

reed_solomon* reed_solomon_new(int data_shards, int parity_shards) {
	  gf* vm = NULL;
    gf* top = NULL;
    int err = 0;
    reed_solomon* rs = NULL;

    assert(fec_initialized);

    do {
        rs = RS_MALLOC(sizeof(reed_solomon));
        if(NULL == rs) {
            return NULL;
        }
        rs->data_shards = data_shards;
        rs->parity_shards = parity_shards;
        rs->shards = (data_shards + parity_shards);
        rs->m = NULL;
        rs->parity = NULL;

        if(rs->shards > DATA_SHARDS_MAX || data_shards <= 0 || parity_shards <= 0) {
            err = 1;
            break;
        }

        vm = vandermonde(rs->shards, rs->data_shards);
        if(NULL == vm) {
            err = 2;
            break;
        }

        top = sub_matrix(vm, 0, 0, data_shards, data_shards, rs->shards, data_shards);
        if(NULL == top) {
            err = 3;
            break;
        }

        err = invert_mat(top, data_shards);
        assert(0 == err);

        rs->m = multiply1(vm, rs->shards, data_shards, top, data_shards, data_shards);
        if(NULL == rs->m) {
            err = 4;
            break;
        }

        rs->parity = sub_matrix(rs->m, data_shards, 0, rs->shards, data_shards, rs->shards, data_shards);
        if(NULL == rs->parity) {
            err = 5;
            break;
        }

        RS_FREE(vm);
        RS_FREE(top);
        vm = NULL;
        top = NULL;
        return rs;

    } while(0);

    fprintf(stderr, "err=%d\n", err);
    if(NULL != vm) {
        RS_FREE(vm);
    }
    if(NULL != top) {
        RS_FREE(top);
    }
    if(NULL != rs) {
        if(NULL != rs->m) {
            RS_FREE(rs->m);
        }
        if(NULL != rs->parity) {
            RS_FREE(rs->parity);
        }
        RS_FREE(rs);
    }

    return NULL;
}

void reed_solomon_release(reed_solomon* rs) {
    if(NULL != rs) {
        if(NULL != rs->m) {
            RS_FREE(rs->m);
        }
        if(NULL != rs->parity) {
            RS_FREE(rs->parity);
        }
        RS_FREE(rs);
    }
}

int reed_solomon_encode(reed_solomon* rs,
        unsigned char** data_blocks,
        unsigned char** fec_blocks,
        int block_size) {
    assert(NULL != rs && NULL != rs->parity);

    return code_some_shards(rs->parity, data_blocks, fec_blocks
            , rs->data_shards, rs->parity_shards, block_size);
}

int reed_solomon_decode(reed_solomon* rs,
        unsigned char **data_blocks,
        int block_size,
        unsigned char **dec_fec_blocks,
        unsigned int *fec_block_nos,
        unsigned int *erased_blocks,
        int nr_fec_blocks) {

    gf dataDecodeMatrix[DATA_SHARDS_MAX*DATA_SHARDS_MAX];
    unsigned char* subShards[DATA_SHARDS_MAX];
    unsigned char* outputs[DATA_SHARDS_MAX];
    gf* m = rs->m;
    int i, j, c, swap, subMatrixRow, dataShards, nos, nshards;

    for(i = 0; i < nr_fec_blocks; i++) {
        swap = 0;
        for(j = i+1; j < nr_fec_blocks; j++) {
            if(erased_blocks[i] > erased_blocks[j]) {

                c = erased_blocks[i];
                erased_blocks[i] = erased_blocks[j];
                erased_blocks[j] = c;

                swap = 1;
            }
        }
        if(!swap) {
            break;
        }
    }

    j = 0;
    subMatrixRow = 0;
    nos = 0;
    nshards = 0;
    dataShards = rs->data_shards;
    for(i = 0; i < dataShards; i++) {
        if(j < nr_fec_blocks && i == erased_blocks[j]) {
            j++;
        } else {
            for(c = 0; c < dataShards; c++) {
                dataDecodeMatrix[subMatrixRow*dataShards + c] = m[i*dataShards + c];
            }
            subShards[subMatrixRow] = data_blocks[i];
            subMatrixRow++;
        }
    }

    for(i = 0; i < nr_fec_blocks && subMatrixRow < dataShards; i++) {
        subShards[subMatrixRow] = dec_fec_blocks[i];
        j = dataShards + fec_block_nos[i];
        for(c = 0; c < dataShards; c++) {
            dataDecodeMatrix[subMatrixRow*dataShards + c] = m[j*dataShards + c];
        }
        subMatrixRow++;
    }

    if(subMatrixRow < dataShards) {
        return -1;
    }

    invert_mat(dataDecodeMatrix, dataShards);

    for(i = 0; i < nr_fec_blocks; i++) {
        j = erased_blocks[i];
        outputs[i] = data_blocks[j];

        memmove(dataDecodeMatrix+i*dataShards, dataDecodeMatrix+j*dataShards, dataShards);
    }

    return code_some_shards(dataDecodeMatrix, subShards, outputs,
            dataShards, nr_fec_blocks, block_size);
}

int reed_solomon_encode2(reed_solomon* rs, unsigned char** shards, int nr_shards, int block_size) {
    unsigned char** data_blocks;
    unsigned char** fec_blocks;
    int i, ds = rs->data_shards, ps = rs->parity_shards, ss = rs->shards;
    i = nr_shards / ss;
    data_blocks = shards;
    fec_blocks = &shards[(i*ds)];

    for(i = 0; i < nr_shards; i += ss) {
        reed_solomon_encode(rs, data_blocks, fec_blocks, block_size);
        data_blocks += ds;
        fec_blocks += ps;
    }
    return 0;
}

int reed_solomon_reconstruct(reed_solomon* rs,
        unsigned char** shards,
        unsigned char* marks,
        int nr_shards,
        int block_size) {
    unsigned char *dec_fec_blocks[DATA_SHARDS_MAX];
    unsigned int fec_block_nos[DATA_SHARDS_MAX];
    unsigned int erased_blocks[DATA_SHARDS_MAX];
    unsigned char* fec_marks;
    unsigned char **data_blocks, **fec_blocks;
    int i, j, dn, pn, n;
    int ds = rs->data_shards;
    int ps = rs->parity_shards;
    int err = 0;

    data_blocks = shards;
    n = nr_shards / rs->shards;
    fec_marks = marks + n*ds;
    fec_blocks = shards + n*ds;

    for(j = 0; j < n; j++) {
        dn = 0;
        for(i = 0; i < ds; i++) {
            if(marks[i]) {
                erased_blocks[dn++] = i;
            }
        }
        if(dn > 0) {
            pn = 0;
            for(i = 0; i < ps && pn < dn; i++) {
                if(!fec_marks[i]) {
                    fec_block_nos[pn] = i;
                    dec_fec_blocks[pn] = fec_blocks[i];
                    pn++;
                }
            }

            if(dn == pn) {
                reed_solomon_decode(rs
                        , data_blocks
                        , block_size
                        , dec_fec_blocks
                        , fec_block_nos
                        , erased_blocks
                        , dn);
            } else {
                err = -1;
            }
        }
        data_blocks += ds;
        marks += ds;
        fec_blocks += ps;
        fec_marks += ps;
    }

    return err;
	}

	info* info_new(char *fname, int dataShards, int parityShards, int blockSize) {
		info* information = malloc( sizeof(info) ) ;
		struct stat st ;
		int i, n, nrFecBlocks ;

		int fd = open(fname, O_RDONLY) ;

		if(fd < 0) {
			fprintf(stderr, "input file : %s not found\n", fname) ;
			exit(1) ;
		}

		information->blockSize = blockSize ;
		information->dataShards = dataShards ;
		information->parityShards = parityShards ;

		fstat(fd, &st) ;
		information->size = st.st_size ;
		information->nrBlocks = (information->size + information->blockSize - 1) / information->blockSize ;
		information->nrBlocks = ((information->nrBlocks + information->dataShards - 1) / information->dataShards) * information->dataShards ;
		i = information->nrBlocks / information->dataShards ;
		information->nrFecBlocks = i * information->parityShards ;
		information->nrShards = information->nrBlocks + information->nrFecBlocks ;

		information->data = (unsigned char*)malloc(information->nrShards * information->blockSize) ;

		n = read(fd, information->data, information->size) ;

		if(n < information->size)
		{
					fprintf(stderr, "Short read\n") ;
					close(fd) ;
					exit(1) ;
		}
		else close(fd) ;

		memset(information->data + information->size, 0, information->nrShards * information->blockSize - information->size) ;

		return information ;
	}

unsigned char** xorENC(reed_solomon* rs, info* inf) {
		unsigned char **data_blocks ;

		data_blocks = (unsigned char**)malloc( inf->nrShards * sizeof(unsigned char**) ) ;
	  for(int i = 0 ; i < inf->nrShards ; i++) {
			data_blocks[i] = inf->data + i * inf->blockSize ;
	  }

		clock_t start, end ;

		start = clock() ;

		reed_solomon_encode2(rs, data_blocks, inf->nrBlocks + inf->nrFecBlocks, inf->blockSize) ;

		end = clock() ;

		fprintf(stderr, "reed_solomon_encode2 execution time: %f\n", (double)(end - start)/CLOCKS_PER_SEC) ;

		return data_blocks ;
}

unsigned char* sim(reed_solomon* rs, info* inf) {
		unsigned char *zilch ;
		int  corrupted ;

		zilch = (unsigned char*)calloc(1, inf->nrShards) ;
	  memset(zilch, 0, sizeof(zilch)) ;
	  corrupted = inf->parityShards ;

		for(int i = 0 ; i < corrupted ; i++) {
        int corr = random() % inf->nrBlocks ;
        memset(inf->data + inf->blockSize * corr, 'E', inf->blockSize) ;
        fprintf(stderr, "Corrupting %d\n", corr) ;
        zilch[corr] = 1 ;
		}

		return zilch ;
}

void xorDEC(reed_solomon* rs, info* inf, unsigned char* zilch, unsigned char** data_blocks) {
		clock_t start, end ;

		start = clock() ;

		reed_solomon_reconstruct(rs, data_blocks, zilch, inf->nrShards, inf->blockSize);

		end = clock() ;

		fprintf(stderr, "reed_solomon_reconstruct execution time: %f\n", (double)(end - start)/CLOCKS_PER_SEC) ;
}
*/
import "C"

import (
	"flag"
	"fmt"
	"os"
	"time"
)

var dataShards = flag.Int("d", 10, "Number of shards to split the data into, must be below 257")
var parityShards = flag.Int("p", 3, "Number of parity shards")
var blockSize = flag.Int("b", 1024, "Size of block")

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
	if *dataShards > 257 {
		fmt.Fprintf(os.Stderr, "Error: Too many data shards\n")
		os.Exit(1)
	}
	fname := args[0]

	C.fec_init()

	startTime := time.Now()

	var rs *C.reed_solomon
	rs = C.reed_solomon_new(C.int(*dataShards), C.int(*parityShards))

	elapsedTime := time.Since(startTime)

	fmt.Fprintf(os.Stderr, "new excution time: %f\n", elapsedTime.Seconds())

	var info *C.info
	info = C.info_new(C.CString(fname), C.int(*dataShards), C.int(*parityShards), C.int(*blockSize))

	fmt.Fprintf(os.Stderr, " size : %d\n blockSize : %d\n nrShards : %d\n nrBlocks : %d\n nrFecBlocks : %d\n dataShards : %d\n parityShards : %d\n", info.size, info.blockSize, info.nrShards, info.nrBlocks, info.nrFecBlocks, info.dataShards, info.parityShards)

	var dataBlocks **C.uchar
	dataBlocks = C.xorENC(rs, info)

	var loss *C.uchar
	loss = C.sim(rs, info)

	C.xorDEC(rs, info, loss, dataBlocks)

	// C.write(1, unsafe.Pointer(info.data), C.ulong(info.size))
}

func checkErr(err error) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err.Error())
		os.Exit(2)
	}
}
