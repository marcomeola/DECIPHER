/****************************************************************************
 *                  Nucleotide compression/decompression                    *
 *                           Author: Erik Wright                            *
 ****************************************************************************/

/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

/* for Calloc/Free */
#include <R_ext/RS.h>

// DECIPHER header file
#include "DECIPHER.h"

// stdint library
#include "stdint.h"

////////////////////////////////////////////////////
// modified CRC calculators from https://pycrc.org/
////////////////////////////////////////////////////

/*
 *    Width      = 8
 *    Poly       = 0x07
 *    XorIn      = 0x00
 *    ReflectIn  = False
 *    XorOut     = 0x00
 *    ReflectOut = False
 */

typedef uint_fast8_t crc_t8;

static const crc_t8 crc_table8[256] = {
	0x00, 0x07, 0x0e, 0x09, 0x1c, 0x1b, 0x12, 0x15, 0x38, 0x3f, 0x36, 0x31, 0x24, 0x23, 0x2a, 0x2d,
	0x70, 0x77, 0x7e, 0x79, 0x6c, 0x6b, 0x62, 0x65, 0x48, 0x4f, 0x46, 0x41, 0x54, 0x53, 0x5a, 0x5d,
	0xe0, 0xe7, 0xee, 0xe9, 0xfc, 0xfb, 0xf2, 0xf5, 0xd8, 0xdf, 0xd6, 0xd1, 0xc4, 0xc3, 0xca, 0xcd,
	0x90, 0x97, 0x9e, 0x99, 0x8c, 0x8b, 0x82, 0x85, 0xa8, 0xaf, 0xa6, 0xa1, 0xb4, 0xb3, 0xba, 0xbd,
	0xc7, 0xc0, 0xc9, 0xce, 0xdb, 0xdc, 0xd5, 0xd2, 0xff, 0xf8, 0xf1, 0xf6, 0xe3, 0xe4, 0xed, 0xea,
	0xb7, 0xb0, 0xb9, 0xbe, 0xab, 0xac, 0xa5, 0xa2, 0x8f, 0x88, 0x81, 0x86, 0x93, 0x94, 0x9d, 0x9a,
	0x27, 0x20, 0x29, 0x2e, 0x3b, 0x3c, 0x35, 0x32, 0x1f, 0x18, 0x11, 0x16, 0x03, 0x04, 0x0d, 0x0a,
	0x57, 0x50, 0x59, 0x5e, 0x4b, 0x4c, 0x45, 0x42, 0x6f, 0x68, 0x61, 0x66, 0x73, 0x74, 0x7d, 0x7a,
	0x89, 0x8e, 0x87, 0x80, 0x95, 0x92, 0x9b, 0x9c, 0xb1, 0xb6, 0xbf, 0xb8, 0xad, 0xaa, 0xa3, 0xa4,
	0xf9, 0xfe, 0xf7, 0xf0, 0xe5, 0xe2, 0xeb, 0xec, 0xc1, 0xc6, 0xcf, 0xc8, 0xdd, 0xda, 0xd3, 0xd4,
	0x69, 0x6e, 0x67, 0x60, 0x75, 0x72, 0x7b, 0x7c, 0x51, 0x56, 0x5f, 0x58, 0x4d, 0x4a, 0x43, 0x44,
	0x19, 0x1e, 0x17, 0x10, 0x05, 0x02, 0x0b, 0x0c, 0x21, 0x26, 0x2f, 0x28, 0x3d, 0x3a, 0x33, 0x34,
	0x4e, 0x49, 0x40, 0x47, 0x52, 0x55, 0x5c, 0x5b, 0x76, 0x71, 0x78, 0x7f, 0x6a, 0x6d, 0x64, 0x63,
	0x3e, 0x39, 0x30, 0x37, 0x22, 0x25, 0x2c, 0x2b, 0x06, 0x01, 0x08, 0x0f, 0x1a, 0x1d, 0x14, 0x13,
	0xae, 0xa9, 0xa0, 0xa7, 0xb2, 0xb5, 0xbc, 0xbb, 0x96, 0x91, 0x98, 0x9f, 0x8a, 0x8d, 0x84, 0x83,
	0xde, 0xd9, 0xd0, 0xd7, 0xc2, 0xc5, 0xcc, 0xcb, 0xe6, 0xe1, 0xe8, 0xef, 0xfa, 0xfd, 0xf4, 0xf3
};

crc_t8 crc_update8(crc_t8 crc, const void *data, int data_len)
{
	const unsigned char *d = (const unsigned char *)data;
	unsigned int tbl_idx;
	
	while (data_len--) {
		tbl_idx = (crc ^ *d);
		crc = (crc_table8[tbl_idx]) & 0xff;
		
		d++;
	}
	
	return crc & 0xff;
}

/*
 *    Width      = 16
 *    Poly       = 0x8005
 *    XorIn      = 0x0000
 *    ReflectIn  = False
 *    XorOut     = 0x0000
 *    ReflectOut = False
 */

typedef uint_fast16_t crc_t16;

static const crc_t16 crc_table16[256] = {
	0x0000, 0xc0c1, 0xc181, 0x0140, 0xc301, 0x03c0, 0x0280, 0xc241,
	0xc601, 0x06c0, 0x0780, 0xc741, 0x0500, 0xc5c1, 0xc481, 0x0440,
	0xcc01, 0x0cc0, 0x0d80, 0xcd41, 0x0f00, 0xcfc1, 0xce81, 0x0e40,
	0x0a00, 0xcac1, 0xcb81, 0x0b40, 0xc901, 0x09c0, 0x0880, 0xc841,
	0xd801, 0x18c0, 0x1980, 0xd941, 0x1b00, 0xdbc1, 0xda81, 0x1a40,
	0x1e00, 0xdec1, 0xdf81, 0x1f40, 0xdd01, 0x1dc0, 0x1c80, 0xdc41,
	0x1400, 0xd4c1, 0xd581, 0x1540, 0xd701, 0x17c0, 0x1680, 0xd641,
	0xd201, 0x12c0, 0x1380, 0xd341, 0x1100, 0xd1c1, 0xd081, 0x1040,
	0xf001, 0x30c0, 0x3180, 0xf141, 0x3300, 0xf3c1, 0xf281, 0x3240,
	0x3600, 0xf6c1, 0xf781, 0x3740, 0xf501, 0x35c0, 0x3480, 0xf441,
	0x3c00, 0xfcc1, 0xfd81, 0x3d40, 0xff01, 0x3fc0, 0x3e80, 0xfe41,
	0xfa01, 0x3ac0, 0x3b80, 0xfb41, 0x3900, 0xf9c1, 0xf881, 0x3840,
	0x2800, 0xe8c1, 0xe981, 0x2940, 0xeb01, 0x2bc0, 0x2a80, 0xea41,
	0xee01, 0x2ec0, 0x2f80, 0xef41, 0x2d00, 0xedc1, 0xec81, 0x2c40,
	0xe401, 0x24c0, 0x2580, 0xe541, 0x2700, 0xe7c1, 0xe681, 0x2640,
	0x2200, 0xe2c1, 0xe381, 0x2340, 0xe101, 0x21c0, 0x2080, 0xe041,
	0xa001, 0x60c0, 0x6180, 0xa141, 0x6300, 0xa3c1, 0xa281, 0x6240,
	0x6600, 0xa6c1, 0xa781, 0x6740, 0xa501, 0x65c0, 0x6480, 0xa441,
	0x6c00, 0xacc1, 0xad81, 0x6d40, 0xaf01, 0x6fc0, 0x6e80, 0xae41,
	0xaa01, 0x6ac0, 0x6b80, 0xab41, 0x6900, 0xa9c1, 0xa881, 0x6840,
	0x7800, 0xb8c1, 0xb981, 0x7940, 0xbb01, 0x7bc0, 0x7a80, 0xba41,
	0xbe01, 0x7ec0, 0x7f80, 0xbf41, 0x7d00, 0xbdc1, 0xbc81, 0x7c40,
	0xb401, 0x74c0, 0x7580, 0xb541, 0x7700, 0xb7c1, 0xb681, 0x7640,
	0x7200, 0xb2c1, 0xb381, 0x7340, 0xb101, 0x71c0, 0x7080, 0xb041,
	0x5000, 0x90c1, 0x9181, 0x5140, 0x9301, 0x53c0, 0x5280, 0x9241,
	0x9601, 0x56c0, 0x5780, 0x9741, 0x5500, 0x95c1, 0x9481, 0x5440,
	0x9c01, 0x5cc0, 0x5d80, 0x9d41, 0x5f00, 0x9fc1, 0x9e81, 0x5e40,
	0x5a00, 0x9ac1, 0x9b81, 0x5b40, 0x9901, 0x59c0, 0x5880, 0x9841,
	0x8801, 0x48c0, 0x4980, 0x8941, 0x4b00, 0x8bc1, 0x8a81, 0x4a40,
	0x4e00, 0x8ec1, 0x8f81, 0x4f40, 0x8d01, 0x4dc0, 0x4c80, 0x8c41,
	0x4400, 0x84c1, 0x8581, 0x4540, 0x8701, 0x47c0, 0x4680, 0x8641,
	0x8201, 0x42c0, 0x4380, 0x8341, 0x4100, 0x81c1, 0x8081, 0x4040
};

crc_t16 crc_update16(crc_t16 crc, const void *data, int data_len)
{
	const unsigned char *d = (const unsigned char *)data;
	unsigned int tbl_idx;
	
	while (data_len--) {
		tbl_idx = (crc ^ *d) & 0xff;
		crc = (crc_table16[tbl_idx] ^ (crc >> 8)) & 0xffff;
		
		d++;
	}
	
	return crc & 0xffff;
}

/*
 *    Width      = 24
 *    Poly       = 0x864cfb
 *    XorIn      = 0xb704ce
 *    ReflectIn  = False
 *    XorOut     = 0x000000
 *    ReflectOut = False
 */

typedef uint_fast32_t crc_t;

static const crc_t crc_table24[256] = {
	0x000000, 0x864cfb, 0x8ad50d, 0x0c99f6, 0x93e6e1, 0x15aa1a, 0x1933ec, 0x9f7f17,
	0xa18139, 0x27cdc2, 0x2b5434, 0xad18cf, 0x3267d8, 0xb42b23, 0xb8b2d5, 0x3efe2e,
	0xc54e89, 0x430272, 0x4f9b84, 0xc9d77f, 0x56a868, 0xd0e493, 0xdc7d65, 0x5a319e,
	0x64cfb0, 0xe2834b, 0xee1abd, 0x685646, 0xf72951, 0x7165aa, 0x7dfc5c, 0xfbb0a7,
	0x0cd1e9, 0x8a9d12, 0x8604e4, 0x00481f, 0x9f3708, 0x197bf3, 0x15e205, 0x93aefe,
	0xad50d0, 0x2b1c2b, 0x2785dd, 0xa1c926, 0x3eb631, 0xb8faca, 0xb4633c, 0x322fc7,
	0xc99f60, 0x4fd39b, 0x434a6d, 0xc50696, 0x5a7981, 0xdc357a, 0xd0ac8c, 0x56e077,
	0x681e59, 0xee52a2, 0xe2cb54, 0x6487af, 0xfbf8b8, 0x7db443, 0x712db5, 0xf7614e,
	0x19a3d2, 0x9fef29, 0x9376df, 0x153a24, 0x8a4533, 0x0c09c8, 0x00903e, 0x86dcc5,
	0xb822eb, 0x3e6e10, 0x32f7e6, 0xb4bb1d, 0x2bc40a, 0xad88f1, 0xa11107, 0x275dfc,
	0xdced5b, 0x5aa1a0, 0x563856, 0xd074ad, 0x4f0bba, 0xc94741, 0xc5deb7, 0x43924c,
	0x7d6c62, 0xfb2099, 0xf7b96f, 0x71f594, 0xee8a83, 0x68c678, 0x645f8e, 0xe21375,
	0x15723b, 0x933ec0, 0x9fa736, 0x19ebcd, 0x8694da, 0x00d821, 0x0c41d7, 0x8a0d2c,
	0xb4f302, 0x32bff9, 0x3e260f, 0xb86af4, 0x2715e3, 0xa15918, 0xadc0ee, 0x2b8c15,
	0xd03cb2, 0x567049, 0x5ae9bf, 0xdca544, 0x43da53, 0xc596a8, 0xc90f5e, 0x4f43a5,
	0x71bd8b, 0xf7f170, 0xfb6886, 0x7d247d, 0xe25b6a, 0x641791, 0x688e67, 0xeec29c,
	0x3347a4, 0xb50b5f, 0xb992a9, 0x3fde52, 0xa0a145, 0x26edbe, 0x2a7448, 0xac38b3,
	0x92c69d, 0x148a66, 0x181390, 0x9e5f6b, 0x01207c, 0x876c87, 0x8bf571, 0x0db98a,
	0xf6092d, 0x7045d6, 0x7cdc20, 0xfa90db, 0x65efcc, 0xe3a337, 0xef3ac1, 0x69763a,
	0x578814, 0xd1c4ef, 0xdd5d19, 0x5b11e2, 0xc46ef5, 0x42220e, 0x4ebbf8, 0xc8f703,
	0x3f964d, 0xb9dab6, 0xb54340, 0x330fbb, 0xac70ac, 0x2a3c57, 0x26a5a1, 0xa0e95a,
	0x9e1774, 0x185b8f, 0x14c279, 0x928e82, 0x0df195, 0x8bbd6e, 0x872498, 0x016863,
	0xfad8c4, 0x7c943f, 0x700dc9, 0xf64132, 0x693e25, 0xef72de, 0xe3eb28, 0x65a7d3,
	0x5b59fd, 0xdd1506, 0xd18cf0, 0x57c00b, 0xc8bf1c, 0x4ef3e7, 0x426a11, 0xc426ea,
	0x2ae476, 0xaca88d, 0xa0317b, 0x267d80, 0xb90297, 0x3f4e6c, 0x33d79a, 0xb59b61,
	0x8b654f, 0x0d29b4, 0x01b042, 0x87fcb9, 0x1883ae, 0x9ecf55, 0x9256a3, 0x141a58,
	0xefaaff, 0x69e604, 0x657ff2, 0xe33309, 0x7c4c1e, 0xfa00e5, 0xf69913, 0x70d5e8,
	0x4e2bc6, 0xc8673d, 0xc4fecb, 0x42b230, 0xddcd27, 0x5b81dc, 0x57182a, 0xd154d1,
	0x26359f, 0xa07964, 0xace092, 0x2aac69, 0xb5d37e, 0x339f85, 0x3f0673, 0xb94a88,
	0x87b4a6, 0x01f85d, 0x0d61ab, 0x8b2d50, 0x145247, 0x921ebc, 0x9e874a, 0x18cbb1,
	0xe37b16, 0x6537ed, 0x69ae1b, 0xefe2e0, 0x709df7, 0xf6d10c, 0xfa48fa, 0x7c0401,
	0x42fa2f, 0xc4b6d4, 0xc82f22, 0x4e63d9, 0xd11cce, 0x575035, 0x5bc9c3, 0xdd8538
};

crc_t crc_update24(crc_t crc, const void *data, int data_len)
{
	const unsigned char *d = (const unsigned char *)data;
	unsigned int tbl_idx;
	
	while (data_len--) {
		tbl_idx = ((crc >> 16) ^ *d) & 0xff;
		crc = (crc_table24[tbl_idx] ^ (crc << 8)) & 0xffffff;
		
		d++;
	}
	
	return crc & 0xffffff;
}

/*
 *    Width      = 32
 *    Poly       = 0x04c11db7
 *    XorIn      = 0xffffffff
 *    ReflectIn  = False
 *    XorOut     = 0x00000000
 *    ReflectOut = False
*/

static const crc_t crc_table32[256] = {
	0x00000000, 0x77073096, 0xee0e612c, 0x990951ba,
	0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3,
	0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
	0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91,
	0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de,
	0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
	0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec,
	0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5,
	0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
	0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,
	0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940,
	0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
	0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116,
	0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f,
	0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
	0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,
	0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a,
	0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
	0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818,
	0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01,
	0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
	0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457,
	0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c,
	0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
	0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2,
	0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb,
	0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
	0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9,
	0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086,
	0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
	0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4,
	0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad,
	0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a,
	0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683,
	0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8,
	0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
	0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe,
	0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7,
	0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
	0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
	0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252,
	0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
	0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60,
	0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79,
	0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
	0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f,
	0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04,
	0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
	0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a,
	0x9c0906a9, 0xeb0e363f, 0x72076785, 0x05005713,
	0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38,
	0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21,
	0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e,
	0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
	0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c,
	0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45,
	0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
	0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db,
	0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0,
	0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
	0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6,
	0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf,
	0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
	0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d
};

crc_t crc_update32(crc_t crc, const void *data, int data_len)
{
	const unsigned char *d = (const unsigned char *)data;
	unsigned int tbl_idx;
	
	while (data_len--) {
		tbl_idx = (crc ^ *d) & 0xff;
		crc = (crc_table32[tbl_idx] ^ (crc >> 8)) & 0xffffffff;
		
		d++;
	}
	
	return crc & 0xffffffff;
}

////////////////////////////////////////////////////
// n-bit compression and decompression algorithms
////////////////////////////////////////////////////

static const unsigned char table[] = {
	0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f,
	0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0x0f,
	0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b,
	0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0x0b,
	0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27,
	0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x07,
	0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23,
	0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x03,
	0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e,
	0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0x0e,
	0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a,
	0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0x0a,
	0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26,
	0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x06,
	0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22,
	0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x02,
	0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d,
	0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0x0d,
	0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29,
	0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x09,
	0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25,
	0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x05,
	0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21,
	0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x01,
	0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c,
	0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0x0c,
	0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28,
	0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x08,
	0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24,
	0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x04,
	0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20,
	0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x00
};

unsigned char revcomp(unsigned char x)
{
	return table[x];
}

unsigned int revcomp2(unsigned int x)
{
	unsigned int temp = 0;
	temp |= table[x & 0xFF] << 8;
	temp |= table[x >> 8];
	return temp;
}

unsigned char reorder(unsigned char byte)
{
	unsigned char x = 0;
	x |= (byte >> 6) & 0x03;
	x |= (byte >> 2) & 0x0C;
	x |= (byte << 2) & 0x30;
	x |= (byte << 6) & 0xC0;
	
	return x;
}

int revcompDiff(const char curr, const char last)
{
	switch (curr) {
		case 'A':
			if (last!='T')
				return 1;
			break;
		case 'C':
			if (last!='G')
				return 1;
			break;
		case 'G':
			if (last!='C')
				return 1;
			break;
		case 'T':
			if (last!='A')
				return 1;
			break;
		case 'U':
			if (last!='A')
				return 1;
			break;
		case 'a':
			if (last!='t')
				return 1;
			break;
		case 'c':
			if (last!='g')
				return 1;
			break;
		case 'g':
			if (last!='c')
				return 1;
			break;
		case 't':
			if (last!='a')
				return 1;
			break;
		case 'u':
			if (last!='a')
				return 1;
			break;
		default:
			return 1;
			break;
	}
	return 0;
}

////////////////////////////////////////////////////
// nbit compression encoding description
//
// First byte:
// 1aab cdee
// aa = 10 for nbit
// b = 0 (start upper case)
// b = 1 (start lower case)
// c = 0 (plain ASCII)
// c = 1 (encoded)
// d = 0 (DNA letters)
// d = 1 (RNA letters)
// ee = number of bytes to store length then CRC:
// 00 = 1 byte (up to 255 nts)
// 01 = 2 bytes (up to 65,535 nts)
// 10 = 3 bytes (up to 16,777,215 nts)
// 11 = 4 bytes (up to 4,294,967,295 nts)
//
// The next bytes provide the length in typical integer format,
// followed by the CRC-8/16/24/32 in accordance with the length.
//
// The next bytes store the sequence in quasi-reverse order:
// A = 00
// C = 01
// G = 10
// T = 00
//
// Examples:
// AAAC = 01000000
// AACA = 00010000
// ACAA = 00000100
// CAAA = 00000001
//
// With the exception of AAAA (0x0), which is the control code.
// When 0x0 is encountered, switch paths based on the next bit:
// axxx xxxx
// a = 0 for path #1
// a = 1 for path #2
//
// Path #1 ('a' = 0):
// If this byte = 0x0 then assign to "AAAA" and continue;
// Otherwise, determine the character from the last 5 bits:
// 0yy0 0001 = A, (length - 1) stored in next byte
// 0yy0 0010 = C, (length - 1) stored in next byte
// 0yy0 0011 = G, (length - 1) stored in next byte
// 0yy0 0100 = T/U, (length - 1) stored in next byte
// 0yy0 0101 = V, (length - 1) stored in next byte
// 0yy0 0110 = H, (length - 1) stored in next byte
// 0yy0 0111 = D, (length - 1) stored in next byte
// 0yy0 1000 = B, (length - 1) stored in next byte
// 0yy0 1001 = +, (length - 1) stored in next byte
// 0yy0 1010 = ., (length - 1) stored in next byte
// 0yy0 1011 = N, (length - 1) stored in next byte
// 0yy0 1100 = -, (length - 1) stored in next byte
// 0yy0 1101 = M, (length - 1) stored in next byte
// 0yy0 1110 = R, (length - 1) stored in next byte
// 0yy0 1111 = W, (length - 1) stored in next byte
// 0yy1 0000 = S, (length - 1) stored in next byte
// 0yy1 0001 = Y, (length - 1) stored in next byte
// 0yy1 0010 = K, (length - 1) stored in next byte
// 0yy1 0011 = a single N
// 0yy1 0100 = a single -
// 0yy1 0101 = a single M
// 0yy1 0110 = a single R
// 0yy1 0111 = a single W
// 0yy1 1000 = a single S
// 0yy1 1001 = a single Y
// 0yy1 1010 = a single K
// 0yy1 1011 = +, (length - 1) stored in next two bytes
// 0yy1 1100 = ., (length - 1) stored in next two bytes
// 0yy1 1101 = N, (length - 1) stored in next two bytes
// 0yy1 1110 = -, (length - 1) stored in next two bytes
// 0yy1 1111 = switch case (upper to lower, or vise-versa)
//
// And the yy bits decide which position to start from:
// yy = 00 for the current position
// yy = 01 for the (current - 1) position
// yy = 10 for the (current - 2) position
// yy = 11 for the (current - 3) position
//
// Path #2 ('a' = 1):
// Three bit encoding until a zero is encountered in the 'a' bit:
// ACTG- = 1-125 combinations (7 bits required)
// 0x0 is still used as the control code (does not break out of triplet coding)
// If the 8-bit is zero then the next byte is no longer triplet code (unless 0x0)
// A = 0, C = 1, G = 2, T = 3, - = 4
// first * 25 + second * 5 + third + 1 = 7-bits
// Examples:
// AAA = 0*25 + 0*5 + 0 + 1 = 1
// --- = 4*25 + 4*5 + 4 + 1 = 125
// A-T = 0*25 + 4*5 + 3 + 1 = 24
//
// 126 & 127 are used as control codes meaning the next 2*n bytes contain:
// start position of the a repeat (n bytes), end position of a repeat
// 126 means that the repeat is exact, 127 means that it is reverse complement
// (Note that the first digit is a 1, so 126 is actually 254 and 127 is 255)
// n = 1 byte, up to position 255
// n = 2 bytes, up to position 65,535
// n = 3 bytes, up to position 16,777,215
// n = 4 bytes, up to position 4,294,967,295
////////////////////////////////////////////////////

// nbit compression algorithm
SEXP nbit(SEXP x, SEXP y, SEXP compRepeats, SEXP nThreads)
{
	int i, j, k, pos;
	int n = length(x);
	int ascii = asInteger(y);
	int cR = asInteger(compRepeats);
	int nthreads = asInteger(nThreads);
	
	unsigned char *p;
	unsigned char **ptrs = Calloc(n, unsigned char *); // compressed strings
	const char *s;
	const char **strs = Calloc(n, const char *); // uncompressed strings
	int *l = Calloc(n, int); // lengths
	
	// build a vector of thread-safe pointers
	for (i = 0; i < n; i++) {
		strs[i] = CHAR(STRING_ELT(x, i));
		l[i] = length(STRING_ELT(x, i));
	}
	
	// compress the sequences
	#pragma omp parallel for private(i,j,k,p,s,pos) schedule(guided) num_threads(nthreads)
	for (i = 0; i < n; i++) {
		ptrs[i] = Calloc(l[i] > 3 ? l[i] : 4, unsigned char); // initialized to zero
		p = ptrs[i];
		s = strs[i];
		
		// initialize the dictionary
		unsigned int *dict, word, count, lastHit, currHit, lastPos = 0;
		int lastTemp, currTemp, rev, len, len2, thresh = 1;
		if (cR==1) {
			if (l[i] <= 5120) { // 20*2^8 = 5120 (~5% chance of an 8-mer occurring only once)
				dict = Calloc(256, unsigned int);
			} else {
				dict = Calloc(65536, unsigned int);
			}
		}
		
		// set 4/7/8-bits to 1
		p[0] |= 200; // 11001000
		
		// set 5-bit (leave as 0 if upper)
		int lower = 0;
		for (j = 0; j < l[i]; j++) {
			if (s[j] >= 'A' && s[j] <= 'Z') {
				break; // uppercase
			} else if (s[j] >= 'a' && s[j] <= 'z') {
				lower = 1; // lowercase
				p[0] |= 16;
				break;
			}
		}
		
		// set 3-bit (leave as 0 if includes T)
		int DNA = 1;
		for (j = 0; j < l[i]; j++) {
			if (s[j]=='T' || s[j]=='t') {
				break; // assume DNA
			} else if (s[j]=='U' || s[j]=='u') {
				DNA = 0;
				p[0] |= 4; // assume RNA
				break;
			}
		}
		
		// set the header length
		int c; // byte count
		// set 1-bit to 2-bit
		if (l[i] > 16777215) {
			p[0] |= 3;
			c = 9;
			p[4] = (l[i] >> 24) & 0xFF;
			p[3] = (l[i] >> 16) & 0xFF;
			p[2] = (l[i] >> 8) & 0xFF;
			p[1] = l[i] & 0xFF;
			
			crc_t crc = 0xffffffff;
			crc = crc_update32(crc, s, l[i]);
			p[8] = (unsigned char)((crc >> 24) & 0xFF);
			p[7] = (unsigned char)((crc >> 16) & 0xFF);
			p[6] = (unsigned char)((crc >> 8) & 0xFF);
			p[5] = (unsigned char)(crc & 0xFF);
		} else if (l[i] > 65535) {
			p[0] |= 2;
			c = 7;
			p[3] = (l[i] >> 16) & 0xFF;
			p[2] = (l[i] >> 8) & 0xFF;
			p[1] = l[i] & 0xFF;
			
			crc_t crc = 0xb704ce;
			crc = crc_update24(crc, s, l[i]);
			p[6] = (unsigned char)((crc >> 16) & 0xFF);
			p[5] = (unsigned char)((crc >> 8) & 0xFF);
			p[4] = (unsigned char)(crc & 0xFF);
		} else if (l[i] > 255) {
			p[0] |= 1;
			c = 5;
			p[2] = (l[i] >> 8) & 0xFF;
			p[1] = l[i] & 0xFF;
			
			crc_t16 crc = 0x0000;
			crc = crc_update16(crc, s, l[i]);
			p[4] = (unsigned char)((crc >> 8) & 0xFF);
			p[3] = (unsigned char)(crc & 0xFF);
		} else {
			c = 3;
			p[1] = l[i] & 0xFF;
			
			crc_t8 crc = 0x00;
			crc = crc_update8(crc, s, l[i]);
			p[2] = (unsigned char)crc;
		}
		
		j = 0; // position in sequence
		pos = 0; // start of run
		int run, lastTriplet, lastCase;
		int threeBit = 0; // 3-bit encoding
		int threeBitBegin = 0; // starting position of 3-bit encoding
		int threeBitEnd = 0; // ending position of 3-bit encoding
		int lastGap = -1; // position of the last gap
		int pair = 0; // pair of bits
		unsigned char byte = 0; // encoded byte
		int success = 1; // successful compression
		while (j < l[i]) {
			run = 0;
			if (j >= pos) { // check run length
				for (k = 1, pos = j + 1; (pos < l[i]) && (k < 256); k++, pos++) {
					if (s[j] != s[pos])
						break;
				}
				if (k > 12) {
					if (s[j]=='+' ||
						s[j]=='.' ||
						s[j]=='-' ||
						s[j]=='N' ||
						s[j]=='n') {
						// look for extended run
						for (pos = j + k; (pos < l[i]) && (k < 65536); k++, pos++) {
							if (s[j] != s[pos])
								break;
						}
					}
					if (pair==0) {
						run = 1;
					} else {
						if (threeBit==0) {
							if (pair==2) {
								pos = j + 3;
							} else if (pair==4) {
								pos = j + 2;
							} else { // pair==6
								pos = j; // check next iteration
							}
						} else {
							if (pair==1) {
								pos = j + 2;
							} else { // pair==2
								pos = j; // check next iteration
							}
						}
					}
				} else if (threeBit==0 && k==1 && s[j]=='-') {
					// look ahead for the next gap
					run = 1; // isolated gap
					for (int h = j + 1; (h < l[i]) && (h <= j + 20); h++) {
						if (s[h]=='-') {
							run = 0; // start 3-bit encoding
							break;
						}
					}
				}
			}
			
			if (s[j]=='A' && run==0) {
				if (lower==1)
					goto switchCase;
				if (threeBit==1) {
					pair++;
				} else {
					pair += 2;
				}
			} else if (s[j]=='a' && run==0) {
				if (lower==0)
					goto switchCase;
				if (threeBit==1) {
					pair++;
				} else {
					pair += 2;
				}
			} else if (s[j]=='C' && run==0) {
				if (lower==1)
					goto switchCase;
				if (threeBit==1) {
					byte += (pair==0) ? 25 : ((pair==1) ? 5:1);
					pair++;
				} else {
					byte |= 1 << pair;
					pair += 2;
				}
			} else if (s[j]=='c' && run==0) {
				if (lower==0)
					goto switchCase;
				if (threeBit==1) {
					byte += (pair==0) ? 25 : ((pair==1) ? 5:1);
					pair++;
				} else {
					byte |= 1 << pair;
					pair += 2;
				}
			} else if (s[j]=='G' && run==0) {
				if (lower==1)
					goto switchCase;
				if (threeBit==1) {
					byte += (pair==0) ? 50 : ((pair==1) ? 10:2);
					pair++;
				} else {
					byte |= 2 << pair;
					pair += 2;
				}
			} else if (s[j]=='g' && run==0) {
				if (lower==0)
					goto switchCase;
				if (threeBit==1) {
					byte += (pair==0) ? 50 : ((pair==1) ? 10:2);
					pair++;
				} else {
					byte |= 2 << pair;
					pair += 2;
				}
			} else if (s[j]=='T' && run==0) {
				if (lower==1)
					goto switchCase;
				if (DNA==0) {
					success = 0;
					break;
				}
				if (threeBit==1) {
					byte += (pair==0) ? 75 : ((pair==1) ? 15:3);
					pair++;
				} else {
					byte |= 3 << pair;
					pair += 2;
				}
			} else if (s[j]=='t' && run==0) {
				if (lower==0)
					goto switchCase;
				if (DNA==0) {
					success = 0;
					break;
				}
				if (threeBit==1) {
					byte += (pair==0) ? 75 : ((pair==1) ? 15:3);
					pair++;
				} else {
					byte |= 3 << pair;
					pair += 2;
				}
			} else if (s[j]=='U' && run==0) {
				if (lower==1)
					goto switchCase;
				if (DNA==1) {
					success = 0;
					break;
				}
				if (threeBit==1) {
					byte += (pair==0) ? 75 : ((pair==1) ? 15:3);
					pair++;
				} else {
					byte |= 3 << pair;
					pair += 2;
				}
			} else if (s[j]=='u' && run==0) {
				if (lower==0)
					goto switchCase;
				if (DNA==1) {
					success = 0;
					break;
				}
				if (threeBit==1) {
					byte += (pair==0) ? 75 : ((pair==1) ? 15:3);
					pair++;
				} else {
					byte |= 3 << pair;
					pair += 2;
				}
			} else if (s[j]=='-' && run==0) {
				if (threeBit==0) {
					threeBit = 1;
					if ((c + 1) >= l[i]) {
						success = 0; // compression failed
						break;
					}
					byte = 129; // 10000001
					p[c] = 0;
					c++;
					threeBitBegin = c;
					if (pair!=0) { // current position
						// unable to reverse in 3-bit encoding
						// retract and switch coding
						j -= (pair==2) ? 1 : ((pair==4) ? 2 : 3);
						pair = 0;
						continue;
					}
				} else if (c > threeBitBegin) {
					lastGap = c;
				}
				
				byte += (pair==0) ? 100 : ((pair==1) ? 20:4);
				pair++;
			} else { // use control code (0)
				// special character or run
				// 3 bytes (nul char/code reps)
				// bit 8 is 0 for runs
				if (s[j] >= 'A' && s[j] <= 'Z') {
					if (lower==1)
						goto switchCase;
				} else if (s[j] >= 'a' && s[j] <= 'z') {
					if (lower==0)
						goto switchCase;
				}
				
				if (pair==0) { // current position
					if ((c + 3) >= l[i]) {
						success = 0; // compression failed
						break;
					}
					p[c] = 0;
					// pair==0 (start from current position)
					byte = 0;
					c++;
				} else { // first record previous positions
					if ((c + 4) >= l[i]) {
						success = 0; // compression failed
						break;
					}
					if (byte==0) // force non-zero byte
						byte |= 1 << 6; // AAAC
					p[c++] = byte;
					p[c++] = 0;
					
					// record starting position
					if (threeBit==0) {
						if (pair==2) { // current - 3
							pair = 96;
						} else if (pair==4) { // current - 2
							pair = 64;
						} else if (pair==6) { // current - 1
							pair = 32;
						}
					} else {
						if (pair==1) { // current - 2
							pair = 64;
						} else if (pair==2) { // current - 1
							pair = 32;
						}
					}
					
					byte = 0;
					byte |= pair; // pair = bits 6/7
				}
				
				int letter;
				switch (s[j]) {
					case 'A':
					case 'a':
						letter = 1;
						break;
					case 'C':
					case 'c':
						letter = 2;
						break;
					case 'G':
					case 'g':
						letter = 3;
						break;
					case 'T':
					case 't':
						if (DNA==0) {
							success = 0;
							break;
						}
						letter = 4;
						break;
					case 'U':
					case 'u':
						if (DNA==1) {
							success = 0;
							break;
						}
						letter = 4;
						break;
					case 'M':
					case 'm':
						letter = 13;
						break;
					case 'R':
					case 'r':
						letter = 14;
						break;
					case 'W':
					case 'w':
						letter = 15;
						break;
					case 'S':
					case 's':
						letter = 16;
						break;
					case 'Y':
					case 'y':
						letter = 17;
						break;
					case 'K':
					case 'k':
						letter = 18;
						break;
					case 'V':
					case 'v':
						letter = 5;
						break;
					case 'H':
					case 'h':
						letter = 6;
						break;
					case 'D':
					case 'd':
						letter = 7;
						break;
					case 'B':
					case 'b':
						letter = 8;
						break;
					case 'N':
					case 'n':
						letter = 11;
						break;
					case '-':
						letter = 12;
						break;
					case '+':
						letter = 9;
						break;
					case '.':
						letter = 10;
						break;
					default:
						success = 0; // compression failed
						break;
				}
				
				if (success==0)
					break;
				
				if (k==1 && letter >= 11) {
					letter += 8;
					byte |= letter;
					p[c] = byte;
					c++;
				} else if (k > 256) {
					letter += 18;
					byte |= letter;
					p[c] = byte;
					c++;
					p[c] = ((k - 1) >> 8) & 0xFF; // length of run
					c++;
					p[c] = (k - 1) & 0xFF; // length of run
					c++;
				} else {
					byte |= letter;
					p[c] = byte;
					c++;
					p[c] = k - 1; // length of run
					c++;
				}
				
				if (threeBit==0) {
					byte = 0;
				} else {
					byte = 129; // 10000001
				}
				j += k;
				pair = 0;
				continue;
			}
			
			if ((pair==8 || j==(l[i] - 1)) && threeBit==0) {
				len = 0;
				if (cR==1) {
					// find previous occurrences of a large region
					// record non-overlapping k-mers in dict
					// look for extendable k-mers in dict
					if (lastPos != (j - 4)) {
						// initialize
						word = (unsigned int)reorder(byte);
						count = 1;
					} else {
						word = (word << 8) | (unsigned int)reorder(byte);
						count++;
						
						// determine the min length required
						if (j > 16777215) {
							thresh = 40; // 10 bytes
						} else if (j > 65535) {
							thresh = 32; // 8 bytes
						} else if (j > 255) {
							thresh = 24; // 6 bytes
						} else {
							thresh = 16; // 4 bytes
						}
						
						if (l[i] <= 5120 && count==2) { // use single byte indices
							// look for repeats in dictionary
							for (k = 0; k < 8; k += 2) {
								currHit = j - 3 - (k >> 1); // start of byte
								
								// exact repeats
								lastHit = dict[(word >> k) & 0xFF];
								if (lastHit!=0) {
									lastHit -= 3; // start of lastHit
									// extend hit
									len = 4;
									for (lastTemp = lastHit + len, currTemp = currHit + len; currTemp < l[i]; len++, lastTemp++, currTemp++) {
										if (s[lastTemp]!=s[currTemp])
											break;
									}
									len -= k >> 1;
									
									if (len >= thresh) {
										// check that 4-mer's case matches
										for (len2 = 1, lastTemp = lastHit + 1; len2 < 4; len2++, lastTemp++) {
											if (s[currHit + len2]!=s[lastTemp]) {
												len = 0;
												break;
											}
										}
										if (len >= thresh) {
											lastHit += k >> 1;
											rev = 0;
											break;
										}
									}
								}
								
								// revcomp repeats
								lastHit = dict[revcomp((word >> k) & 0xFF)]; // end of lastHit
								if (lastHit!=0) {
									// extend hit
									len = 4;
									for (lastTemp = lastHit - len, currTemp = currHit + len; lastTemp >= 0; len++, lastTemp--, currTemp++) {
										if (revcompDiff(s[currTemp], s[lastTemp]))
											break;
									}
									len -= k >> 1;
									
									if (len >= thresh) {
										// check that 4-mer's case matches
										for (len2 = 1, lastTemp = lastHit - 1; len2 < 4; len2++, lastTemp--) {
											if (revcompDiff(s[currHit + len2], s[lastTemp])) {
												len = 0;
												break;
											}
										}
										if (len >= thresh) {
											lastHit -= k >> 1;
											rev = 1;
											break;
										}
									}
								}
							}
							
							// record starting position in dictionary
							word = word & 0xFF;
							dict[word] = j;
							count = 1;
						} else if (l[i] > 5120 && count==4) { // use double byte indices
							// look for repeats in dictionary
							for (k = 0; k < 16; k += 2) {
								currHit = j - 3 - (k >> 1); // start of byte
								
								// exact repeats
								lastHit = dict[(word >> k) & 0xFFFF];
								if (lastHit!=0) {
									lastHit -= 3; // start of lastHit
									// extend hit
									len = 4;
									for (lastTemp = lastHit + len, currTemp = currHit + len; currTemp < l[i]; len++, lastTemp++, currTemp++) {
										if (s[lastTemp]!=s[currTemp])
											break;
									}
									len -= k >> 1;
									
									if (len >= thresh) {
										// check that 4-mer's case matches
										for (len2 = 1, lastTemp = lastHit + 1; len2 < 4; len2++, lastTemp++) {
											if (s[currHit + len2]!=s[lastTemp]) {
												len = 0;
												break;
											}
										}
										if (len >= thresh) {
											lastHit += k >> 1;
											rev = 0;
											break;
										}
									}
								}
								
								// revcomp repeats
								lastHit = dict[revcomp2((word >> k) & 0xFFFF)]; // end of lastHit
								if (lastHit!=0) {
									// extend hit
									len = 4;
									for (lastTemp = lastHit - len - 4, currTemp = currHit + len; lastTemp >= 0; len++, lastTemp--, currTemp++) {
										if (revcompDiff(s[currTemp], s[lastTemp]))
											break;
									}
									len -= k >> 1;
									
									if (len >= thresh) {
										// check that 4-mer's case matches
										for (len2 = 1, lastTemp = lastHit - 5; len2 < 4; len2++, lastTemp--) {
											if (revcompDiff(s[currHit + len2], s[lastTemp])) {
												len = 0;
												break;
											}
										}
										if (len >= thresh) {
											lastHit -= (k >> 1) + 4;
											rev = 1;
											break;
										}
									}
								}
							}
							
							// record starting position in dictionary
							word = word & 0xFFFF;
							dict[word & 0xFFFF] = j;
							count = 2;
						}
					}
					lastPos = j;
				}
				
				if (len >= thresh) { // repeat
					j -= 3; // start of byte
					
					if (j > 16777215) {
						if ((c + 9) >= l[i]) {
							success = 0; // compression failed
							break;
						}
						p[c++] = 0;
						p[c++] = rev==0 ? 254 : 255;
						p[c++] = (unsigned char)(lastHit >> 24);
						p[c++] = (unsigned char)(lastHit >> 16);
						p[c++] = (unsigned char)(lastHit >> 8);
						p[c++] = (unsigned char)lastHit;
						p[c++] = (unsigned char)(len >> 24);
						p[c++] = (unsigned char)(len >> 16);
						p[c++] = (unsigned char)(len >> 8);
						p[c++] = (unsigned char)len;
					} else if (j > 65535) {
						if ((c + 7) >= l[i]) {
							success = 0; // compression failed
							break;
						}
						if (len > 16777215)
							len = 16777215;
						p[c++] = 0;
						p[c++] = rev==0 ? 254 : 255;
						p[c++] = (unsigned char)(lastHit >> 16);
						p[c++] = (unsigned char)(lastHit >> 8);
						p[c++] = (unsigned char)lastHit;
						p[c++] = (unsigned char)(len >> 16);
						p[c++] = (unsigned char)(len >> 8);
						p[c++] = (unsigned char)len;
					} else if (j > 255) {
						if ((c + 5) >= l[i]) {
							success = 0; // compression failed
							break;
						}
						if (len > 65535)
							len = 65535;
						p[c++] = 0;
						p[c++] = rev==0 ? 254 : 255;
						p[c++] = (unsigned char)(lastHit >> 8);
						p[c++] = (unsigned char)lastHit;
						p[c++] = (unsigned char)(len >> 8);
						p[c++] = (unsigned char)len;
					} else {
						if ((c + 3) >= l[i]) {
							success = 0; // compression failed
							break;
						}
						if (len > 255)
							len = 255;
						p[c++] = 0;
						p[c++] = rev==0 ? 254 : 255;
						p[c++] = (unsigned char)lastHit;
						p[c++] = (unsigned char)len;
					}
					
					j += len - 1;
					byte = 0;
					pair = 0;
				} else if (byte==0) { // AAAA
					if ((c + 1) >= l[i]) {
						success = 0; // compression failed
						break;
					}
					// repeat zero byte twice
					p[c++] = 0;
					p[c++] = 0;
				} else {
					if (c >= l[i] && c > 2) {
						success = 0; // compression failed
						break;
					}
					p[c] = byte;
					c++;
					byte = 0;
				}
				pair = 0;
			} else if (pair==3 || j==(l[i] - 1)) {
				if (threeBitEnd > threeBitBegin && (j - lastTriplet) > 20) {
					// re-encode using 2-bit encoding because it is
					// more efficient (20/3 ~= 20/4 + 1 + partial byte)
					p[threeBitEnd] &= 127; // zero 8-bit
					c = threeBitEnd + 1;
					j = lastTriplet + 1;
					pos = j;
					threeBit = 0;
					byte = 0;
					pair = 0;
					lower = lastCase;
					continue;
				}
				if (c >= l[i] && c > 2) {
					success = 0; // compression failed
					break;
				}
				if (c==lastGap) {
					threeBitEnd = c;
					lastTriplet = j;
					lastCase = lower;
				}
				
				p[c] = byte;
				c++;
				byte = 129; // 10000001
				pair=0;
			}
			
			j++;
			continue;
			
			switchCase:
			if (lower==0) {
				lower = 1;
			} else {
				lower = 0;
			}
			
			if (pair==0) {
				if ((c + 2) >= l[i]) {
					success = 0; // compression failed
					break;
				}
				p[c++] = 0;
				p[c] = 31;
				// byte is already correct
			} else {
				if ((c + 3) >= l[i]) {
					success = 0; // compression failed
					break;
				}
				
				if (threeBit==0) {
					if (byte==0) // force non-zero byte
						byte |= 1 << 6; // AAAC
					p[c++] = byte;
					byte = 0;
					p[c++] = 0;
					
					if (pair==2) { // current - 3
						p[c] = 127;
					} else if (pair==4) { // current - 2
						p[c] = 95;
					} else if (pair==6) { // current - 1
						p[c] = 63;
					}
				} else {
					p[c++] = byte;
					byte = 129; // 10000001
					p[c++] = 0;
					
					if (pair==1) { // current - 2
						p[c] = 95;
					} else if (pair==2) { // current - 1
						p[c] = 63;
					}
				}
				pair = 0;
			}
			
			pos = j;
			c++;
		}
		
		if (cR==1)
			Free(dict);
		
		if (success==0) {
			l[i] = 0;
			p[0] &= 247; // zero the 4-bit
		} else {
			l[i] = c; // new length
		}
	}
	
	Free(strs);
	
	SEXP ret, ans;
	PROTECT(ret = allocVector(VECSXP, n));
	
	for (i = 0; i < n; i++) {
		p = ptrs[i];
		if (l[i]==0) { // compression failed
			if (ascii==1) { // keep as ascii
				l[i] = length(STRING_ELT(x, i));
				PROTECT(ans = allocVector(RAWSXP, l[i] + 1));
				// copy header byte
				RAW(ans)[0] = p[0];
				// copy characters directly
				memcpy(RAW(ans) + 1, CHAR(STRING_ELT(x, i)), l[i]);
			} else { // return empty raw vector
				PROTECT(ans = allocVector(RAWSXP, l[i]));
			}
		} else { // compression succeeded
			PROTECT(ans = allocVector(RAWSXP, l[i]));
			memcpy(RAW(ans), p, l[i]);
		}
		
		Free(p);
		SET_VECTOR_ELT(ret, i, ans);
		UNPROTECT(1); // ans
	}
	
	Free(ptrs);
	Free(l);
	UNPROTECT(1); // ret
	
	return ret;
}

////////////////////////////////////////////////////
// qbit compression encoding description
//
// First byte:
// 1aab bbcc
// aa = 01 for qbit
// bbb = minimum difference between levels (delta)
// cc = number of bytes to store length then CRC:
// 00 = 1 byte (up to 255 nts)
// 01 = 2 bytes (up to 65,535 nts)
// 10 = 3 bytes (up to 16,777,215 nts)
// 11 = 4 bytes (up to 4,294,967,295 nts)
//
// The next bytes provide the length in typical integer format,
// followed by the CRC-8/16/24/32 in accordance with the length.
//
// The following seven bits store the starting value.
//
// The remaining bits store the transformed quality scores:
// (1) Record the value of the first element
// (2) Take the lagged difference of order one
// (3) Apply a bijection mapping to positive integers
// (4) Find the minimum value (offset) greater than one
// (5) Subtract (offset - 1) from all values greater than zero
// (6) Convert runs of >= 8 ones into 0xFF followed by length - 8
// (7) Convert longer runs into 0xFFFF followed by length - 8
// (8) Encode the remainder using Elias gamma encoding
////////////////////////////////////////////////////

// qbit compression algorithm
SEXP qbit(SEXP x, SEXP y, SEXP nThreads)
{
	int i, j;
	int n = length(x);
	int ascii = asInteger(y);
	int nthreads = asInteger(nThreads);
	
	unsigned char leading1[9] = {0, 128, 192, 224, 240, 248, 252, 254, 255};
	unsigned char trailing1[9] = {255, 127, 63, 31, 15, 7, 3, 1, 0};
	unsigned char zeros[256] = {
		0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
	};
	
	unsigned char *p;
	unsigned char **ptrs = Calloc(n, unsigned char *); // compressed strings
	const char *s;
	const char **strs = Calloc(n, const char *); // uncompressed strings
	int *l = Calloc(n, int); // lengths
	
	// build a vector of thread-safe pointers
	for (i = 0; i < n; i++) {
		strs[i] = CHAR(STRING_ELT(x, i));
		l[i] = length(STRING_ELT(x, i));
	}
	
	// compress the quality scores
	#pragma omp parallel for private(i,j,k,p,s,pos) schedule(guided) num_threads(nthreads)
	for (i = 0; i < n; i++) {
		ptrs[i] = Calloc(l[i] > 3 ? l[i] : 4, unsigned char); // initialized to zero
		p = ptrs[i];
		s = strs[i];
		
		int success = 1; // successful compression
		for (j = 0; j < l[i]; j++) {
			if (s[j] > 126) {
				success = 0;
				break;
			}
		}
		
		// set the header length
		int c; // byte count
		// set 1-bit to 2-bit
		if (l[i] > 16777215) {
			p[0] |= 3;
			c = 9;
			p[4] = (l[i] >> 24) & 0xFF;
			p[3] = (l[i] >> 16) & 0xFF;
			p[2] = (l[i] >> 8) & 0xFF;
			p[1] = l[i] & 0xFF;
			
			crc_t crc = 0xffffffff;
			crc = crc_update32(crc, s, l[i]);
			p[8] = (unsigned char)((crc >> 24) & 0xFF);
			p[7] = (unsigned char)((crc >> 16) & 0xFF);
			p[6] = (unsigned char)((crc >> 8) & 0xFF);
			p[5] = (unsigned char)(crc & 0xFF);
		} else if (l[i] > 65535) {
			p[0] |= 2;
			c = 7;
			p[3] = (l[i] >> 16) & 0xFF;
			p[2] = (l[i] >> 8) & 0xFF;
			p[1] = l[i] & 0xFF;
			
			crc_t crc = 0xb704ce;
			crc = crc_update24(crc, s, l[i]);
			p[6] = (unsigned char)((crc >> 16) & 0xFF);
			p[5] = (unsigned char)((crc >> 8) & 0xFF);
			p[4] = (unsigned char)(crc & 0xFF);
		} else if (l[i] > 255) {
			p[0] |= 1;
			c = 5;
			p[2] = (l[i] >> 8) & 0xFF;
			p[1] = l[i] & 0xFF;
			
			crc_t16 crc = 0x0000;
			crc = crc_update16(crc, s, l[i]);
			p[4] = (unsigned char)((crc >> 8) & 0xFF);
			p[3] = (unsigned char)(crc & 0xFF);
		} else {
			c = 3;
			p[1] = l[i] & 0xFF;
			
			crc_t8 crc = 0x00;
			crc = crc_update8(crc, s, l[i]);
			p[2] = (unsigned char)crc;
		}
		
		if (success) {
			// set 6/8-bits to 1
			p[0] |= 160; // 10100000
		} else {
			// use the nbit header
			l[i] = 0;
			p[0] |= 192; // 11000000
			continue;
		}
		
		// fill the next 7 bits with the first character
		if (l[i] > 0) {
			p[c] = (s[0] << 1) & 0xFF;
		} else {
			continue;
		}
		
		char min = 127; // minimum value > 1
		int temp;
		l[i]--;
		unsigned char *t = Calloc(l[i], unsigned char);
		for (j = 0; j < l[i]; j++) {
			temp = s[j + 1] - s[j];
			if (temp > 0) {
				t[j] = 2*temp;
			} else {
				t[j] = -2*temp + 1;
			}
			if (t[j] > 1 && t[j] < min)
				min = t[j];
		}
		
		// record offset and remove
		if (min != 127 && min != 2) {
			min -= 2;
			if (min > 7)
				min = 7;
			p[0] |= min << 2;
			
			// subtract min from values > 1
			for (j = 0; j < l[i]; j++) {
				if (t[j] > 1)
					t[j] -= min;
			}
		}
		
		int b = 7; // current bit in byte
		j = 0; // position in sequence
		int k = 0; // length of run
		int pos = 0; // start of run
		unsigned char byte;
		while (j < l[i]) {
			if (j >= pos && t[j] == 1) { // check run length
				for (k = 1, pos = j + 1; (pos < l[i]) && (k < 65544); k++, pos++) {
					if (t[j] != t[pos])
						break;
				}
			}
			
			if (k > 7) { // run of ones
				if (k < 263) {
					if ((c + 2) > l[i]) {
						success = 0;
						break;
					}
					
					// add eight ones
					p[c] |= trailing1[b];
					c++;
					p[c] |= leading1[b];
					
					byte = k - 8;
					p[c] |= byte >> b;
					c++;
					p[c] |= byte << (8 - b);
				} else {
					if ((c + 4) > l[i]) {
						success = 0;
						break;
					}
					
					// add eight ones
					p[c] |= trailing1[b];
					c++;
					p[c] = 255;
					c++;
					p[c] |= leading1[b];
					
					byte = (k - 8) >> 8;
					p[c] |= byte >> b;
					c++;
					p[c] |= byte << (8 - b);
					byte = k - 8;
					p[c] |= byte >> b;
					c++;
					p[c] |= byte << (8 - b);
				}
				
				j += k;
				k = 0;
				continue;
			}
			
			// apply Elias gamma encoding
			b += zeros[t[j]];
			if (b > 7) {
				c++;
				if (c > l[i]) {
					success = 0;
					break;
				}
				b -= 8;
			}
			
			temp = 7 - b - zeros[t[j]];
			if (temp == 0) {
				p[c] |= t[j];
				c++;
				if (c > l[i]) {
					success = 0;
					break;
				}
				b = 0;
			} else if (temp > 0) {
				p[c] |= t[j] << temp;
				b += zeros[t[j]] + 1;
			} else {
				temp *= -1;
				p[c] |= t[j] >> temp;
				b += zeros[t[j]] + 1;
				c++;
				if (c > l[i]) {
					success = 0;
					break;
				}
				b -= 8;
				p[c] |= t[j] << (8 - temp);
			}
			
			j++;
		}
		
		Free(t);
		
		if (success==0) {
			l[i] = 0;
			p[0] |= 64; // make the 7-bit one
			p[0] &= 223; // make the 6-bit zero
		} else {
			// set the new length
			if (b==0) {
				l[i] = c;
			} else {
				l[i] = c + 1;
			}
		}
	}
	
	Free(strs);
	
	SEXP ret, ans;
	PROTECT(ret = allocVector(VECSXP, n));
	
	for (i = 0; i < n; i++) {
		p = ptrs[i];
		if (l[i]==0) { // compression failed
			if (ascii==1) { // keep as ascii
				l[i] = length(STRING_ELT(x, i));
				PROTECT(ans = allocVector(RAWSXP, l[i] + 1));
				// copy header byte
				RAW(ans)[0] = p[0];
				// copy characters directly
				memcpy(RAW(ans) + 1, CHAR(STRING_ELT(x, i)), l[i]);
			} else { // return empty raw vector
				PROTECT(ans = allocVector(RAWSXP, l[i]));
			}
		} else { // compression succeeded
			PROTECT(ans = allocVector(RAWSXP, l[i]));
			memcpy(RAW(ans), p, l[i]);
		}
		
		Free(p);
		SET_VECTOR_ELT(ret, i, ans);
		UNPROTECT(1); // ans
	}
	
	Free(ptrs);
	Free(l);
	UNPROTECT(1); // ret
	
	return ret;
}

// decompression algorithm
SEXP decompress(SEXP x, SEXP nThreads)
{
	int i;
	int n = length(x);
	int nthreads = asInteger(nThreads);
	
	unsigned char bit[8] = {128, 64, 32, 16, 8, 4, 2, 1};
	unsigned char bijection[256] = {
		0, 0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9,
		10, -10, 11, -11, 12, -12, 13, -13, 14, -14, 15, -15, 16, -16, 17, -17, 18, -18, 19, -19,
		20, -20, 21, -21, 22, -22, 23, -23, 24, -24, 25, -25, 26, -26, 27, -27, 28, -28, 29, -29,
		30, -30, 31, -31, 32, -32, 33, -33, 34, -34, 35, -35, 36, -36, 37, -37, 38, -38, 39, -39,
		40, -40, 41, -41, 42, -42, 43, -43, 44, -44, 45, -45, 46, -46, 47, -47, 48, -48, 49, -49,
		50, -50, 51, -51, 52, -52, 53, -53, 54, -54, 55, -55, 56, -56, 57, -57, 58, -58, 59, -59,
		60, -60, 61, -61, 62, -62, 63, -63, 64, -64, 65, -65, 66, -66, 67, -67, 68, -68, 69, -69,
		70, -70, 71, -71, 72, -72, 73, -73, 74, -74, 75, -75, 76, -76, 77, -77, 78, -78, 79, -79,
		80, -80, 81, -81, 82, -82, 83, -83, 84, -84, 85, -85, 86, -86, 87, -87, 88, -88, 89, -89,
		90, -90, 91, -91, 92, -92, 93, -93, 94, -94, 95, -95, 96, -96, 97, -97, 98, -98, 99, -99,
		100, -100, 101, -101, 102, -102, 103, -103, 104, -104, 105, -105, 106, -106, 107, -107, 108, -108, 109, -109,
		110, -110, 111, -111, 112, -112, 113, -113, 114, -114, 115, -115, 116, -116, 117, -117, 118, -118, 119, -119,
		120, -120, 121, -121, 122, -122, 123, -123, 124, -124, 125, -125, 126, -126, 127, -127
	};
	
	char *s;
	char **strs = Calloc(n, char *); // uncompressed strings
	unsigned char *p;
	unsigned char **ptrs = Calloc(n, unsigned char *); // compressed strings
	int *l = Calloc(n, int); // lengths
	
	// build a vector of thread-safe pointers
	for (i = 0; i < n; i++) {
		ptrs[i] = RAW(VECTOR_ELT(x, i));
		l[i] = length(VECTOR_ELT(x, i));
		if (l[i]==0)
			error("x contains an empty raw vector.");
	}
	
	#pragma omp parallel for private(i,p,s) schedule(guided) num_threads(nthreads)
	for (i = 0; i < n; i++) {
		p = ptrs[i];
		
		char TU, tu, A, C, G, T;
		int len, j, k, lower;
		unsigned char type, min;
		
		type = (p[0] & 224);
		
		if (type==192 || type==160) { // nbit or qbit compression
			if (type==192 && (p[0] & 8)==0) { // ascii
				s = Calloc(l[i], char); // each sequence
				memcpy(s, p + 1, l[i] - 1);
				strs[i] = s;
				continue;
			} else {
				len = 0;
				if ((p[0] & 3)==3) {
					len |= p[1];
					len |= p[2] << 8;
					len |= p[3] << 16;
					len |= p[4] << 24;
					j = 9;
				} else if ((p[0] & 3)==2) {
					len |= p[1];
					len |= p[2] << 8;
					len |= p[3] << 16;
					j = 7;
				} else if ((p[0] & 3)==1) {
					len |= p[1];
					len |= p[2] << 8;
					j = 5;
				} else {
					len |= p[1];
					j = 3;
				}
			}
			if (type==192) { // nbit compression
				if ((p[0] & 4) > 0) {
					TU = 'U';
					tu = 'u';
				} else {
					TU = 'T';
					tu = 't';
				}
				
				lower = p[0] & 16; // > 0 if lower
				if (lower==0) {
					A = 'A';
					C = 'C';
					G = 'G';
					T = TU;
				} else {
					A = 'a';
					C = 'c';
					G = 'g';
					T = tu;
				}
			} else { // qbit compression
				min = (p[0] & 28) >> 2;
			}
		} else {
			l[i] = 0;
			continue;
		}
		
		// initialized to zero
		s = Calloc(len + 4, char); // each sequence
		strs[i] = s;
		
		// decompress the payload
		if (type==192) { // nbit compression
			int run;
			char letter;
			int byte;
			int c = 0;
			int threeBit = 0;
			while (j < l[i]) {
				if (p[j]==0) { // control code
					j++;
					if (j==l[i])
						error("Corrupted encoding.");
					if (p[j]==0) { // AAAA
						s[c++] = A;
						s[c++] = A;
						s[c++] = A;
						s[c++] = A;
					} else if (p[j]==254 || p[j]==255) { // repeat
						int rev = (p[j]==254) ? 0 : 1;
						unsigned int start = 0;
						unsigned int len = 0;
						j++;
						
						if (c > 16777215) {
							start |= p[j++] << 24;
							start |= p[j++] << 16;
							start |= p[j++] << 8;
							start |= p[j++];
							len |= p[j++] << 24;
							len |= p[j++] << 16;
							len |= p[j++] << 8;
							len |= p[j];
						} else if (c > 65535) {
							start |= p[j++] << 16;
							start |= p[j++] << 8;
							start |= p[j++];
							len |= p[j++] << 16;
							len |= p[j++] << 8;
							len |= p[j];
						} else if (c > 255) {
							start |= p[j++] << 8;
							start |= p[j++];
							len |= p[j++] << 8;
							len |= p[j];
						} else {
							start |= p[j++];
							len |= p[j];
						}
						
						if (rev==0) { // exact repeat
							for (k = 0; k < len; k++, start++, c++)
								s[c] = s[start];
						} else { // revcomp repeat
							for (k = 0; k < len; k++, start--, c++) {
								switch (s[start]) {
									case 'A':
										s[c] = 'T';
										break;
									case 'a':
										s[c] = 't';
										break;
									case 'C':
										s[c] = 'G';
										break;
									case 'c':
										s[c] = 'g';
										break;
									case 'G':
										s[c] = 'C';
										break;
									case 'g':
										s[c] = 'c';
										break;
									case 'T':
										s[c] = 'A';
										break;
									case 't':
										s[c] = 'a';
										break;
									case 'U':
										s[c] = 'A';
										break;
									case 'u':
										s[c] = 'a';
										break;
								}
							}
						}
					} else if (p[j]==31) {
						goto switchCase;
					} else if (p[j]==63) {
						c--;
						goto switchCase;
					} else if (p[j]==95) {
						c -= 2;
						goto switchCase;
					} else if (p[j]==127) {
						c -= 3;
						goto switchCase;
					} else if ((p[j] >> 7)==0) { // run
						c -= p[j] >> 5;
						byte = (p[j] & 31);
						if (byte > 18 && byte < 27) {
							switch (byte) {
								case 21:
									if (lower==0) {
										s[c] = 'M';
									} else {
										s[c] = 'm';
									}
									break;
								case 22:
									if (lower==0) {
										s[c] = 'R';
									} else {
										s[c] = 'r';
									}
									break;
								case 23:
									if (lower==0) {
										s[c] = 'W';
									} else {
										s[c] = 'w';
									}
									break;
								case 24:
									if (lower==0) {
										s[c] = 'S';
									} else {
										s[c] = 's';
									}
									break;
								case 25:
									if (lower==0) {
										s[c] = 'Y';
									} else {
										s[c] = 'y';
									}
									break;
								case 26:
									if (lower==0) {
										s[c] = 'K';
									} else {
										s[c] = 'k';
									}
									break;
								case 19:
									if (lower==0) {
										s[c] = 'N';
									} else {
										s[c] = 'n';
									}
									break;
								case 20:
									s[c] = '-';
									break;
								default:
									error("Unexpected byte.");
									break;
							}
							c++;
						} else if (byte > 26) {
							switch (byte) {
								case 27:
									letter = '+';
									break;
								case 28:
									letter = '.';
									break;
								case 29:
									if (lower==0) {
										letter = 'N';
									} else {
										letter = 'n';
									}
									break;
								case 30:
									letter = '-';
									break;
								default:
									error("Unexpected byte.");
									break;
							}
							j++;
							if (j==l[i])
								error("Corrupted encoding.");
							run = p[j] << 8;
							j++;
							if (j==l[i])
								error("Corrupted encoding.");
							run |= p[j];
							for (k = 0; k <= run; k++, c++) {
								s[c] = letter;
							}
						} else {
							switch (byte) {
								case 1:
									letter = A;
									break;
								case 2:
									letter = C;
									break;
								case 3:
									letter = G;
									break;
								case 4:
									letter = T;
									break;
								case 9:
									letter = '+';
									break;
								case 10:
									letter = '.';
									break;
								case 13:
									if (lower==0) {
										letter = 'M';
									} else {
										letter = 'm';
									}
									break;
								case 14:
									if (lower==0) {
										letter = 'R';
									} else {
										letter = 'r';
									}
									break;
								case 15:
									if (lower==0) {
										letter = 'W';
									} else {
										letter = 'w';
									}
									break;
								case 16:
									if (lower==0) {
										letter = 'S';
									} else {
										letter = 's';
									}
									break;
								case 17:
									if (lower==0) {
										letter = 'Y';
									} else {
										letter = 'y';
									}
									break;
								case 18:
									if (lower==0) {
										letter = 'K';
									} else {
										letter = 'k';
									}
									break;
								case 5:
									if (lower==0) {
										letter = 'V';
									} else {
										letter = 'v';
									}
									break;
								case 6:
									if (lower==0) {
										letter = 'H';
									} else {
										letter = 'h';
									}
									break;
								case 7:
									if (lower==0) {
										letter = 'D';
									} else {
										letter = 'd';
									}
									break;
								case 8:
									if (lower==0) {
										letter = 'B';
									} else {
										letter = 'b';
									}
									break;
								case 11:
									if (lower==0) {
										letter = 'N';
									} else {
										letter = 'n';
									}
									break;
								case 12:
									letter = '-';
									break;
								default:
									error("Unexpected byte.");
									break;
							}
							j++;
							if (j==l[i])
								error("Corrupted encoding.");
							run = p[j];
							for (k = 0; k <= run; k++, c++) {
								s[c] = letter;
							}
						}
					} else { // 3-bit encoding
						threeBit = 1;
						continue; // don't increment j
					}
				} else if (threeBit) { // 3-bit encoding
					if ((p[j] & 128)==0) {
						threeBit = 0; // next byte is not 3-bit encoding
						byte = p[j];
					} else {
						byte = p[j] & 127; // clear 8-bit
					}
					
					if (byte > 100) {
						s[c++] = '-';
						byte -= 100;
					} else if (byte > 75) {
						s[c++] = T;
						byte -= 75;
					} else if (byte > 50) {
						s[c++] = G;
						byte -= 50;
					} else if (byte > 25) {
						s[c++] = C;
						byte -= 25;
					} else {
						s[c++] = A;
					}
					if (byte > 20) {
						s[c++] = '-';
						byte -= 20;
					} else if (byte > 15) {
						s[c++] = T;
						byte -= 15;
					} else if (byte > 10) {
						s[c++] = G;
						byte -= 10;
					} else if (byte > 5) {
						s[c++] = C;
						byte -= 5;
					} else {
						s[c++] = A;
					}
					if (byte==5) {
						s[c++] = '-';
					} else if (byte==4) {
						s[c++] = T;
					} else if (byte==3) {
						s[c++] = G;
					} else if (byte==2) {
						s[c++] = C;
					} else {
						s[c++] = A;
					}
				} else { // 2-bit encoding
					switch (p[j] & 3) {
						case 0:
							s[c++] = A;
							break;
						case 1:
							s[c++] = C;
							break;
						case 2:
							s[c++] = G;
							break;
						case 3:
							s[c++] = T;
							break;
					}
					switch (p[j] & 12) {
						case 0:
							s[c++] = A;
							break;
						case 4:
							s[c++] = C;
							break;
						case 8:
							s[c++] = G;
							break;
						case 12:
							s[c++] = T;
							break;
					}
					switch (p[j] & 48) {
						case 0:
							s[c++] = A;
							break;
						case 16:
							s[c++] = C;
							break;
						case 32:
							s[c++] = G;
							break;
						case 48:
							s[c++] = T;
							break;
					}
					switch (p[j] & 192) {
						case 0:
							s[c++] = A;
							break;
						case 64:
							s[c++] = C;
							break;
						case 128:
							s[c++] = G;
							break;
						case 192:
							s[c++] = T;
							break;
					}
				}
				
				j++;
				continue;
				
				switchCase:
				if (lower==0) {
					lower = 1;
					A = 'a';
					C = 'c';
					G = 'g';
					T = tu;
				} else {
					lower = 0;
					A = 'A';
					C = 'C';
					G = 'G';
					T = TU;
				}
				j++;
			}
		} else { // qbit compression
			s[0] = (p[j] & 254) >> 1;
			
			unsigned char *t = Calloc(len, unsigned char); // t-gaps
			
			int c = 0; // position in t
			int b = 7; // current bit in byte
			int ones = 0;
			int zeros = 0;
			int byte;
			while (j < l[i]) {
				if ((p[j] & bit[b])) { // one
					if (zeros) { // previously zero
						// record value
						zeros++; // number of bits in value
						if ((b + zeros) > 8) { // straddles bytes
							byte = (unsigned char)(p[j] << b) >> (8 - zeros);
							j++;
							if (j==l[i])
								error("Corrupted encoding.");
							b += zeros - 8;
							byte |= p[j] >> (8 - b);
						} else {
							byte = (unsigned char)(p[j] << b) >> (8 - zeros);
							b += zeros;
							if (b==8) {
								b = 0;
								j++;
							}
						}
						t[c++] = byte;
						zeros = 0;
					} else {
						t[c++] = 1;
						if (b==7) {
							b = 0;
							j++;
						} else {
							b++;
						}
						ones++;
						if (ones==8) {
							// run of ones
							if (j==l[i])
								error("Corrupted encoding.");
							byte = (unsigned char)(p[j] << b);
							j++;
							if (b > 0) {
								if (j==l[i])
									error("Corrupted encoding.");
								byte |= p[j] >> (8 - b);
							}
							if (byte==255) {
								byte = (unsigned char)(p[j] << b);
								j++;
								if (j==l[i])
									error("Corrupted encoding.");
								byte |= p[j] >> (8 - b);
								byte <<= 8;
								byte |= (unsigned char)(p[j] << b);
								j++;
								if (b > 0) {
									if (j==l[i])
										error("Corrupted encoding.");
									byte |= p[j] >> (8 - b);
								}
							}
							for (k = 0; k < byte; k++)
								t[c++] = 1;
							ones = 0;
						}
					}
				} else { // zero
					zeros++;
					if (ones)
						ones = 0;
					if (b==7) {
						b = 0;
						j++;
					} else {
						b++;
					}
				}
			}
			
			if (min > 0) {
				for (k = 0; k < c; k++)
					if (t[k] > 1)
						t[k] += min;
			}
			
			// apply the reverse bijection
			for (k = 1; k < len; k++)
				s[k] = s[k - 1] + bijection[t[k - 1]];
			
			Free(t);
		}
		s[len] = '\0'; // null-terminate
		
		// Cyclic Redundancy Check
		if ((p[0] & 3)==3) {
			crc_t crc = 0xffffffff;
			crc = crc_update32(crc, s, len);
			if (p[8] != (unsigned char)((crc >> 24) & 0xFF) ||
				p[7] != (unsigned char)((crc >> 16) & 0xFF) ||
				p[6] != (unsigned char)((crc >> 8) & 0xFF) ||
				p[5] != (unsigned char)(crc & 0xFF))
				l[i] = -1;
		} else if ((p[0] & 3)==2) {
			crc_t crc = 0xb704ce;
			crc = crc_update24(crc, s, len);
			if (p[6] != (unsigned char)((crc >> 16) & 0xFF) ||
				p[5] != (unsigned char)((crc >> 8) & 0xFF) ||
				p[4] != (unsigned char)(crc & 0xFF))
				l[i] = -1;
		} else if ((p[0] & 3)==1) {
			crc_t16 crc = 0x0000;
			crc = crc_update16(crc, s, len);
			if (p[4] != (unsigned char)((crc >> 8) & 0xFF) ||
				p[3] != (unsigned char)(crc & 0xFF))
				l[i] = -1;
		} else {
			crc_t8 crc = 0x00;
			crc = crc_update8(crc, s, len);
			if (p[2] != (unsigned char)crc)
				l[i] = -1;
		}
	}
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	
	for (i = 0; i < n; i++) {
		if (l[i] < 0) {
			error("Data corruption in x[[%d]]", i + 1);
		} else if (l[i]==0) { // not decompressed
			SET_STRING_ELT(seqs, i, NA_STRING);
		} else { // decompressed
			s = strs[i];
			SET_STRING_ELT(seqs, i, mkChar(s));
			Free(s);
		}
	}
	
	Free(ptrs);
	Free(strs);
	Free(l);
	
	UNPROTECT(1);
	
	return seqs;
}
