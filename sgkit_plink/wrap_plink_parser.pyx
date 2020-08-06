import numpy as np 


cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "./CPlinkBedFile.h":

	void _readPlinkBedFilefloatFAAA "readPlinkBedFilefloatFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _readPlinkBedFiledoubleFAAA "readPlinkBedFiledoubleFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _readPlinkBedFilefloatCAAA "readPlinkBedFilefloatCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _readPlinkBedFiledoubleCAAA "readPlinkBedFiledoubleCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _readPlinkBedFileint8FAAA "readPlinkBedFileint8FAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, signed char* out)
	void _readPlinkBedFileint8CAAA "readPlinkBedFileint8CAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, signed char* out)

	void _writePlinkBedFilefloatFAAA "writePlinkBedFilefloatFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, float* inx)
	void _writePlinkBedFiledoubleFAAA "writePlinkBedFiledoubleFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, double* inx)
	void _writePlinkBedFileint8FAAA "writePlinkBedFileint8FAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, signed char* inx)
	void _writePlinkBedFilefloatCAAA "writePlinkBedFilefloatCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, float* inx)
	void _writePlinkBedFiledoubleCAAA "writePlinkBedFiledoubleCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, double* inx)
	void _writePlinkBedFileint8CAAA "writePlinkBedFileint8CAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, signed char* inx)


 
def readPlinkBedFile2floatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=2] out):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <float*> out.data)
	return out

def readPlinkBedFile2floatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=2] out):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <float*> out.data)
	return out


def readPlinkBedFile2doubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=2] out):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <double*> out.data)
	return out

def readPlinkBedFile2doubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=2] out):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <double*> out.data)
	return out

def readPlinkBedFile2int8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.int8_t, ndim=2] out):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFileint8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <signed char*> out.data)
	return out

def readPlinkBedFile2int8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.int8_t, ndim=2] out):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFileint8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <signed char*> out.data)
	return out

def writePlinkBedFile2floatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float32_t, ndim=2] inx):
	_writePlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <float*> inx.data)

def writePlinkBedFile2floatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float32_t, ndim=2] inx):
	_writePlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <float*> inx.data)

def writePlinkBedFile2doubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float64_t, ndim=2] inx):
	_writePlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <double*> inx.data)

def writePlinkBedFile2doubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float64_t, ndim=2] inx):
	_writePlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <double*> inx.data)

def writePlinkBedFile2int8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.int8_t, ndim=2] inx):
	_writePlinkBedFileint8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <signed char*> inx.data)

def writePlinkBedFile2int8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.int8_t, ndim=2] inx):
	_writePlinkBedFileint8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <signed char*> inx.data)

