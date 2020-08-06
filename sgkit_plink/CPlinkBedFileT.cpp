/*
*******************************************************************
*
*    Copyright (c) Microsoft. All rights reserved.
*
*    THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
*    ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
*    IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
*    PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
*
******************************************************************
*/

/*
* CPlinkBedFile - {PLINK BED File Access Class}
*
*         File Name:   CPlinkBedFile.cpp
*           Version:   2.00
*            Author:   
*     Creation Date:   18 Nov 2010
*     Revision Date:   14 Aug 2013
*
*    Module Purpose:   This file implements the CPlinkBedFile 
*  
*                      A .BED file contains compressed binary genotype 
*                         values for individuals by SNPs.  
*
*    Change History:   Version 2.00: Reworked to be wrapped in python version by Chris Widmer  (chris@shogun-toolbox.org)
*
* Test Files: 
*/

/*
* Include Files
*/
#include "CPlinkBedFileT.h"
#include <iostream>
#include <stdio.h>
#include <math.h> 
#include <stdlib.h>

#ifdef _WIN32
#define isinf(x) (!_finite(x))
#else
#define isinf(x) (!isfinite(x))
#endif

#ifdef MISSING_VALUE
REAL SUFFIX(unknownOrMissing) = MISSING_VALUE;
#else
REAL SUFFIX(unknownOrMissing) = std::numeric_limits<REAL>::quiet_NaN();  // now used by SnpInfo
#endif


REAL SUFFIX(homozygousPrimaryAllele) = 0;                // Major Allele
REAL SUFFIX(heterozygousAllele) = 1;                     
REAL SUFFIX(homozygousSecondaryAllele) = 2;              // Minor Allele ()

REAL SUFFIX(mapBedGenotypeToRealAlleleCountA1)[4] = { 
	SUFFIX(homozygousSecondaryAllele),       // look-up 0
	SUFFIX(unknownOrMissing),                // look-up 1
	SUFFIX(heterozygousAllele),              // look-up 2
	SUFFIX(homozygousPrimaryAllele),         // look-up 3
};

REAL SUFFIX(mapBedGenotypeToRealAlleleNoCountA1)[4] = {
	SUFFIX(homozygousPrimaryAllele),         // look-up 0
	SUFFIX(unknownOrMissing),                // look-up 1
	SUFFIX(heterozygousAllele),              // look-up 2
	SUFFIX(homozygousSecondaryAllele),       // look-up 3
};

SUFFIX(CBedFile)::SUFFIX(CBedFile)()
{
	layout = LayoutUnknown;    // layout describes the matrix layout on disk
	// 0=RowMajor(all snps per individual together);
	// 1=ColumnMajor(all individuals per SNP together in memory)
	cIndividuals = 0;
	cSnps        = 0;
	cbStride     = 0;

}

SUFFIX(CBedFile)::SUFFIX(~CBedFile)()
{
	if ( pFile )
	{
		fclose( pFile );
		pFile = NULL;
	}
}


void SUFFIX(CBedFile)::Open( const string& filename_, size_t cIndividuals_, size_t cSnps_ )
{
	if ( filename_.empty() )
	{
		printf( "Could not create BedFile Reader.  Parameter 'filename' is zero length string" );
	}

	filename = filename_;         // TODO: removed FullPath
	cIndividuals = cIndividuals_;
	cSnps = cSnps_;

	pFile = fopen( filename.c_str(), "rb" );  // read in binary to ensure ftell works right
	if ( !pFile )
	{
		printf( "Cannot open input file [%s].\n", filename.c_str()); //TODO: removed errorNO
	}

	//  Verify 'magic' number
	unsigned char rd1 = NextChar();
	unsigned char rd2 = NextChar();
	if ( (bedFileMagic1 != rd1) || (bedFileMagic2 != rd2))
	{
		printf( "Ill-formed BED file [%s]."
			"\n  BED file header is incorrect."
			"\n  Expected magic number of 0x%02x 0x%02x, found 0x%02x 0x%02x", 
			filename.c_str(), bedFileMagic1, bedFileMagic2, rd1, rd2 );
	}

	// Verify 'mode' is valid
	unsigned char rd3 = NextChar();
	switch( rd3 )
	{
	case 0:  // mode = 'IndividualMajor' or RowMajor
		layout = LayoutGroupGenotypesByIndividual;   // all SNPs per individual are sequential in memory
		cbStride = (cSnps + 3)/4;                    // 4 genotypes per byte so round up
		break;
	case 1:  // mode = 'SnpMajor' or ColumnMajor
		layout = LayoutGroupGenotypesBySnp;          // all individuals per SNP are sequential in memory
		cbStride = (cIndividuals + 3)/4;             // 4 genotypes per byte so round up
		break;
	default:
		printf( "Ill-formed BED file [%s].  BED file header is incorrect.  Expected mode to be 0 or 1, found %d", filename.c_str(), rd3 );
		break;
	}

	// allocate the read buffer for a SNP
	rgBytes.resize( cbStride );
	rgBedGenotypes.resize( cIndividuals, bedMissingGenotype );
}

LayoutMode  SUFFIX(CBedFile)::GetLayoutMode()
{
	return( layout );
}

int SUFFIX(CBedFile)::NextChar()
{
	int value = fgetc( pFile );
	if ( value == EOF )
	{
		printf( "Ill-formed BED file [%s]. Encountered EOF before expected.", filename.c_str() );
	}
	return( (unsigned char)value );
}

size_t SUFFIX(CBedFile)::Read( BYTE *pb, size_t cbToRead )
{
	size_t cbRead = fread( pb, 1, cbToRead, pFile );
	if ( cbRead != cbToRead )
	{
		if ( feof( pFile ) )
		{
			printf( "Encountered EOF before expected in BED file. Ill-formed BED file [%s]", filename.c_str() );
		}
		int err = ferror( pFile );
		if ( err )
		{
			printf( "Encountered a file error %d in BED file [%s]", err, filename.c_str() );
		}
	}
	return( cbRead );
}

size_t SUFFIX(CBedFile)::ReadLine(BYTE *pb, size_t idx)
{
	long long fpos = cbHeader + (idx*cbStride);
#ifdef _WIN32
	long long fposCur = _ftelli64(pFile);
#elif __APPLE__
	long long fposCur = ftello(pFile);
#else
	long long fposCur = ftello64(pFile);
#endif
	if (fpos != fposCur)
	{
#ifdef _WIN32
		_fseeki64(pFile, fpos, SEEK_SET);
#elif __APPLE__
		fseeko(pFile, fpos, SEEK_SET);
#else
		fseeko64(pFile, fpos, SEEK_SET);
#endif
	}

	size_t cbRead = Read(pb, cbStride);
	return(cbRead);
}

/*
* Read the genotype for all the individuals in iidList at the SNP specified by iSNP
*   and store the results in pvOut
*/
void SUFFIX(CBedFile)::ReadGenotypes(size_t iSnp, bool count_A1, const vector< size_t >& idxIndividualList, REAL* pvOut, uint64_t_ startpos, uint64_t_  outputNumSNPs)
{
	//fprintf(stdout,"reading iSnp=%d w/ cIndividuals=%d and startpos=%d\n",iSnp,cIndividuals,startpos);
	ReadLine( &rgBytes[0], iSnp );
	// 'decompress' the genotype information
	size_t iIndividual = 0;
	for ( size_t ib = 0; ib < cbStride; ++ib )
	{
		BYTE genotypeByte = rgBytes[ ib ];

		// manually unrolled loop to decompress this byte
		if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)( genotypeByte       & 0x03);
		if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)((genotypeByte >> 2) & 0x03);
		if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)((genotypeByte >> 4) & 0x03);
		if ( iIndividual < cIndividuals ) rgBedGenotypes[ iIndividual++ ] = (BedGenotype)((genotypeByte >> 6) & 0x03);
	}
	for ( size_t i=0; i<idxIndividualList.size(); ++i )
	{
		size_t idx = idxIndividualList[ i ];
		//fprintf(stdout,"iSnp=%d, iIID=%d\n",iSnp,idx);
#ifdef ORDERF
		uint64_t_ out_idx = startpos + i;
#else
		uint64_t_ out_idx = startpos + i * outputNumSNPs;
#endif
		if (count_A1)
		{
			pvOut[out_idx] = SUFFIX(mapBedGenotypeToRealAlleleCountA1)[rgBedGenotypes[idx]];
		}
		else {
			pvOut[out_idx] = SUFFIX(mapBedGenotypeToRealAlleleNoCountA1)[rgBedGenotypes[idx]];
		}

	}
}

// wrapper to be used from cython
void SUFFIX(readPlinkBedFile)(std::string bed_fn, int inputNumIndividuals, int inputNumSNPs, bool count_A1, std::vector<size_t> individuals_idx, std::vector<int> snpIdxList, REAL* out)
{
	uint64_t_ outputNumSNPs = snpIdxList.size();

	SUFFIX(CBedFile) bedFile = SUFFIX(CBedFile)();
	bedFile.Open(bed_fn, inputNumIndividuals, inputNumSNPs);

	for (size_t i = 0; i != snpIdxList.size(); i++) {
		int idx = snpIdxList[i];

#ifdef ORDERF
		uint64_t_ startpos = ((uint64_t_)i) * individuals_idx.size();
#else
		uint64_t_ startpos = ((uint64_t_)i);
#endif
		bedFile.ReadGenotypes(idx, count_A1, individuals_idx, out, startpos, outputNumSNPs);
	}
}

// wrapper to be used from cython
void SUFFIX(writePlinkBedFile)(std::string bed_fn, int iid_count, int sid_count, bool count_A1, REAL* in)
{
	FILE* bed_filepointer = fopen(bed_fn.c_str(), "wb");
	if (!bed_filepointer)
	{
		printf("Cannot open input file [%s].\n", bed_fn.c_str()); //TODO: removed errorNO
		return;
	}

	unsigned char zeroCode = (count_A1 ? 3 : 0);
	unsigned char twoCode  = (count_A1 ? 0 : 3);

	putc(bedFileMagic1, bed_filepointer);
	putc(bedFileMagic2, bed_filepointer);
	putc(1, bed_filepointer);

	//printf("c\n");

	uint64_t_ startpos = 0;
#ifdef ORDERF
	long long int sid_increment = (long long int)0;
	long long int iid_increment = (long long int)1;
#else
	long long int sid_increment = (long long int)1 - (long long int) iid_count*(long long int)sid_count;
	long long int iid_increment = (long long int)sid_count;
#endif

	//printf("d\n");

	for (int sid_index = 0; sid_index < sid_count; ++sid_index)
	{
		//printf("a %d of %d\n", sid_index, sid_count);
		for (int iid_by_four = 0; iid_by_four < iid_count; iid_by_four += 4)
		{
			//printf("e %d of %d\n", iid_by_four, iid_count);
			unsigned char b = 0;

			int end = iid_count - iid_by_four;
			if (end > 4)
			{
				end = 4;
			}
				

			for (int val_index = 0; val_index < end; ++val_index)
			{
				//printf("f %d and %d\n", startpos, sid_index * iid_count + iid_by_four + val_index);
				//printf("g %d of %d\n", val_index, end);
				REAL val = in[startpos];
				unsigned char code;
				if (val == 0)
					code = zeroCode;
				else if (val == 1)
					code = 2; //0b10 backwards on purpose
				else if (val == 2)
					code = twoCode;
#ifdef MISSING_VALUE
				else if (val == MISSING_VALUE)
#else
				else if (val != val)
#endif
					code = 1; //0b01 #backwards on purpose
				else
				{
					//printf("Can't convert value '%s' to BED format (only 0,1,2,NAN[-1] allowed)", val);
					fclose(bed_filepointer);
					return;
				}

				//printf("before b %d, val_index=%d, code=%d, (code << (val_index * 2))=%d\n", b, val_index, code, (code << (val_index * 2)));
				b |= (code << (val_index * 2));
				//printf("code %d makes b %d\n", code, b);
				startpos += iid_increment;
			}
			//printf("writing byte %d\n", b);
			putc(b, bed_filepointer);
		}
		startpos += sid_increment;
		}
	fclose(bed_filepointer);
	//printf("b \n");
}

