#!!!cmk todo: Offer to ignore some or all fam bim fields
#!!!cmk todo: fix up write
#!!!cmk it must be an error to give an override for a property that doesn't exist
#!!!cmk would be nice if the "counts don't match" errors were more specific

#!!!cmk add typing info
#!!!cmk run flake8, isort, etc
import os
import numpy as np
import numbers
import pandas as pd
import logging  #!!!cmk how does sgkit do logging messages?
from pathlib import Path
import multiprocessing

# import warnings
import math
from typing import Any, List, Optional, Tuple, Union

from itertools import takewhile, repeat

# https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
def rawincount(filename):
    f = open(filename, "rb")
    bufgen = takewhile(lambda x: x, (f.raw.read(1024 * 1024) for _ in repeat(None)))
    return sum(buf.count(b"\n") for buf in bufgen)


#!!!cmk test that it works on files with no rows and/or cols (it did before)
class open_bed:  #!!!cmk need doc strings everywhere
    def __init__(
        self,
        filename,
        iid_count=None,
        sid_count=None,
        metadata={},
        count_A1=True,
        num_threads=None,
        skip_format_check=False,
    ):  #!!!document these new optionals. they are here
        self.filename = str(
            Path(filename)
        )  #!!!cmk make it path on the inside, but path or string on the outside
        self.count_A1 = count_A1
        self._num_threads = num_threads
        self.skip_format_check = skip_format_check

        self.metadata_dict, self._counts = open_bed._fixup_metadata(
            metadata, iid_count, sid_count, use_fill_sequence=False
        )
        self._iid_range = None
        self._sid_range = None

        if not self.skip_format_check:
            bedfile = self._name_of_other_file(self.filename, "bed", "bed")
            with open(bedfile, "rb") as filepointer:
                self._check_file(filepointer)

    @staticmethod
    def _all_same(key, length, missing, dtype):
        return np.fill(length, missing, dtype=dtype)

    @staticmethod
    def _sequence(key, length, missing, dtype):
        result = np.zeros(length, dtype=dtype)
        result[:] = (f"{key}{i+1}" for i in range(length))

    _meta_meta = {
        "fid": ("fam", 0, np.str_, "0", _all_same),
        "iid": ("fam", 1, np.str_, None, _sequence),
        "father": ("fam", 2, np.str_, "0", _all_same),
        "mother": ("fam", 3, np.str_, "0", _all_same),
        "sex": ("fam", 4, "int32", 0, _all_same),
        "pheno": ("fam", 5, np.str_, "0", _sequence),
        "chromosome": ("bim", 0, np.str_, "0", _all_same),
        "sid": ("bim", 1, np.str_, None, _sequence),
        "cm_position": ("bim", 2, "float32", 0, _all_same),
        "bp_position": ("bim", 3, "int32", 0, _all_same),
        "allele_1": ("bim", 4, np.str_, "A1", _all_same),
        "allele_2": ("bim", 5, np.str_, "A2", _all_same),
    }

    _delimiters = {"fam": r"\s+", "bim": "\t"}
    _count_name = {"fam": "iid_count", "bim": "sid_count"}

    def _fixup(self, input):
        pass

    @staticmethod
    def _fixup_metadata(metadata, iid_count, sid_count, use_fill_sequence):

        metadata_dict = {key: None for key in open_bed._meta_meta}
        count_dict = {"fam": iid_count, "bim": sid_count}

        for key, input in metadata.items():
            if key not in open_bed._meta_meta:
                raise KeyError(f"metadata key '{key}' not known")

        for (
            key,
            (suffix, _, dtype, missing_value, fill_sequence),
        ) in open_bed._meta_meta.items():
            count = count_dict[suffix]

            input = metadata.get(key)
            if input is None:
                if use_fill_sequence:
                    input = fill_sequence(key, count, missing_value, dtype)
                else:
                    continue
            elif len(input) == 0:
                output = np.zeros([0], dtype=dtype)
            else:
                if not isinstance(input, np.ndarray) or input.dtype.type is not dtype:
                    input = np.array(input, dtype=dtype)
                if len(input.shape) != 1:
                    raise ValueError(f"Override {key} should be one dimensional")
                output = input

            if count is None:
                count_dict[suffix] = len(output)
            else:
                if count != len(output):
                    raise ValueError(
                        f"The length of override {key}, {len(output)}, should not be different from the current {open_bed._count_name[suffix]}, {count}"
                    )
            metadata_dict[key] = output
        return metadata_dict, count_dict

    def _read_fam_or_bim(self, suffix):
        metafile = open_bed._name_of_other_file(self.filename, "bed", suffix)
        logging.info("Loading {0} file {1}".format(suffix, metafile))

        count = self._counts[suffix]

        delimiter = self._delimiters[suffix]
        if delimiter in {r"\s+"}:
            delimiter = None
            delim_whitespace = True
        else:
            delim_whitespace = False

        if os.path.getsize(metafile) == 0:
            fields = []
        else:
            fields = pd.read_csv(
                metafile,
                delimiter=delimiter,
                delim_whitespace=delim_whitespace,
                header=None,
                index_col=False,
                comment=None,
            )

        if count is None:
            self._counts[suffix] = len(fields)
        else:
            if count != len(fields):
                raise ValueError(
                    f"The number of lines in the *.{suffix} file, {len(fields)}, should not be different from the current {self._count_name[suffix]}, {count}"
                )
        #!!!cmk rename suffix_x and other *x variables
        for key, (suffix_x, column, dtype, missing, _) in self._meta_meta.items():
            if suffix_x is not suffix:
                continue
            val = self.metadata_dict[key]
            if val is None:
                if len(fields) == 0:
                    output = np.array([], dtype=dtype)
                if missing is None:
                    output = np.array(fields[column], dtype=dtype)
                else:
                    output = np.array(fields[column].fillna(missing), dtype=dtype)
                self.metadata_dict[key] = output

    @staticmethod
    def _name_of_other_file(filename, remove_suffix, add_suffix):
        if filename.lower().endswith(remove_suffix.lower()):
            filename = filename[0 : -1 - len(remove_suffix)]
        return filename + "." + add_suffix

    def __str__(self):
        return f"{self.__class__.__name__}('{self.filename}',...)"

    @property
    def fid(self):
        return self.metadata_item("fid")

    @property
    def iid(self):
        return self.metadata_item("iid")

    @property
    def father(self):
        return self.metadata_item("father")

    @property
    def mother(self):
        return self.metadata_item("mother")

    @property
    def sex(self):
        return self.metadata_item("sex")

    @property
    def pheno(self):
        return self.metadata_item("pheno")

    @property
    def metadata(self):
        for key in self._meta_meta:
            self.metadata_item(key)
        return self.metadata_dict

    def metadata_item(self, key):
        val = self.metadata_dict.get(key)
        if val is None:
            suffix, _, _, _, _ = self._meta_meta[key]
            self._read_fam_or_bim(suffix=suffix)
            return self.metadata_dict[key]
        else:
            return val

    @property
    def chromosome(self):
        return self.metadata_item("chromosome")

    @property
    def sid(self):
        return self.metadata_item("sid")

    @property
    def cm_position(self):
        return self.metadata_item("cm_position")

    @property
    def bp_position(self):
        return self.metadata_item("bp_position")

    @property
    def allele_1(self):
        return self.metadata_item("allele_1")

    @property
    def allele_2(self):
        return self.metadata_item("allele_2")

    @property
    def iid_count(self):
        return self._count("fam")

    @property
    def sid_count(self):
        return self._count("bim")

    def _count(self, suffix):
        count = self._counts[suffix]
        if count is None:
            metafile = open_bed._name_of_other_file(self.filename, "bed", suffix)
            count = rawincount(metafile)
            self._counts[suffix] = count
        return count

    @staticmethod
    def _check_file(filepointer):
        mode = filepointer.read(2)
        if mode != b"l\x1b":
            raise ValueError(
                "No valid binary BED file"
            )  #!!!cmk make all and any "Exception" more specific
        mode = filepointer.read(1)  # \x01 = SNP major \x00 = individual major
        if mode != b"\x01":
            raise ValueError(
                "only SNP-major is implemented"
            )  #!!!cmk should mention this

    def __del__(self):
        self.__exit__()

    def close(self):
        """
        !!!cmk doc this
            """
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    #!!!cmk say something about support for snp-minor vs major
    @staticmethod
    def write(
        filename, val, metadata, count_A1=True, force_python_only=False,
    ):
        iid_count = val.shape[0]
        sid_count = val.shape[1]

        metadata, _ = open_bed._fixup_metadata(
            metadata, iid_count=iid_count, sid_count=sid_count, use_fill_sequence=False,
        )

        open_bed._write_fam_or_bim(filename, metadata, "fam")
        open_bed._write_fam_or_bim(filename, metadata, "bim")

        bedfile = open_bed._name_of_other_file(
            filename, remove_suffix="bed", add_suffix="bed"
        )

        if not force_python_only:
            from sgkit_plink import wrap_plink_parser

            if val.flags["C_CONTIGUOUS"]:
                order = "C"
            elif val.flags["F_CONTIGUOUS"]:
                order = "F"
            else:
                raise Exception("order not known (not 'F' or 'C')")

            iid_count, sid_count = val.shape
            if val.dtype == np.float64:
                if order == "F":
                    wrap_plink_parser.writePlinkBedFile2doubleFAAA(
                        bedfile.encode("ascii"), iid_count, sid_count, count_A1, val,
                    )
                else:
                    wrap_plink_parser.writePlinkBedFile2doubleCAAA(
                        bedfile.encode("ascii"), iid_count, sid_count, count_A1, val,
                    )
            elif val.dtype == np.float32:
                if order == "F":
                    wrap_plink_parser.writePlinkBedFile2floatFAAA(
                        bedfile.encode("ascii"), iid_count, sid_count, count_A1, val,
                    )
                else:
                    wrap_plink_parser.writePlinkBedFile2floatCAAA(
                        bedfile.encode("ascii"), iid_count, sid_count, count_A1, val,
                    )
            elif val.dtype == np.int8:  #!!!cmk move this up
                if order == "F":
                    wrap_plink_parser.writePlinkBedFile2int8FAAA(
                        bedfile.encode("ascii"), iid_count, sid_count, count_A1, val,
                    )
                else:
                    wrap_plink_parser.writePlinkBedFile2int8CAAA(
                        bedfile.encode("ascii"), iid_count, sid_count, count_A1, val,
                    )
            else:
                raise Exception(
                    "dtype '{0}' not known, only float64 and float32 (and sometimes int8)".format(
                        val.dtype
                    )
                )

        else:
            if not count_A1:
                zero_code = 0b00
                two_code = 0b11
            else:
                zero_code = 0b11
                two_code = 0b00

            with open(bedfile, "wb") as bed_filepointer:
                # see http://zzz.bwh.harvard.edu/plink/binary.shtml
                bed_filepointer.write(bytes(bytearray([0b01101100])))  # magic numbers
                bed_filepointer.write(bytes(bytearray([0b00011011])))  # magic numbers
                bed_filepointer.write(bytes(bytearray([0b00000001])))  # snp major

                for sid_index in range(sid_count):
                    if sid_index % 1 == 0:
                        logging.info(
                            "Writing snp # {0} to file '{1}'".format(
                                sid_index, filename
                            )
                        )

                    colx = val[:, sid_index]
                    for iid_by_four in range(0, len(iid), 4):
                        vals_for_this_byte = colx[iid_by_four : iid_by_four + 4]
                        byte = 0b00000000
                        for val_index in range(len(vals_for_this_byte)):
                            valx = vals_for_this_byte[
                                val_index
                            ]  #!!!cmk rename valx and the other *x
                            if valx == 0:
                                code = zero_code
                            elif valx == 1:
                                code = 0b10  # backwards on purpose
                            elif valx == 2:
                                code = two_code
                            elif (val.dtype == np.int8 and valx == -127) or np.isnan(
                                valx
                            ):  #!!!cmk find a better way to tell int types from float types
                                code = 0b01  # backwards on purpose
                            else:
                                raise ValueError(
                                    "Can't convert value '{0}' to BED format (only 0,1,2,NAN [or sometimes -127] allowed)".format(
                                        valx
                                    )
                                )
                            byte |= code << (val_index * 2)
                        bed_filepointer.write(bytes(bytearray([byte])))
        logging.info("Done writing " + filename)

    def _get_num_threads(self):
        if self._num_threads is not None:
            return self._num_threads
        if "MKL_NUM_THREADS" in os.environ:
            return int(os.environ["MKL_NUM_THREADS"])
        return multiprocessing.cpu_count()

    @staticmethod
    def _array_properties_are_ok(val, order, dtype):
        dtype = np.dtype(dtype)

        if val.dtype != dtype:
            return False
        if order == "F":
            return val.flags["F_CONTIGUOUS"]
        elif order == "C":
            return val.flags["C_CONTIGUOUS"]

        return True

    @property
    def shape(self):
        return (len(self.iid), len(self.sid))

    def read(
        self,
        index: Optional[Any] = None,
        dtype: Optional[Union[type, str]] = np.int8,
        order: Optional[str] = "F",
        force_python_only: bool = False,
    ) -> np.ndarray:

        iid_index_or_none, sid_index_or_none = self._split_index(index)

        if order == "A":
            order = "F"
        dtype = np.dtype(dtype)

        assert not hasattr(
            self, "ind_used"
        ), "A SnpReader should not have a 'ind_used' attribute"

        iid_count_in = (
            self.iid_count
        )  #!!!cmk do we really need all six of these values?
        sid_count_in = self.sid_count

        #!!!cmk happy with _iid_range and _sid_range?
        if self._iid_range is None:
            self._iid_range = np.arange(self.iid_count)
        if self._sid_range is None:
            self._sid_range = np.arange(self.sid_count)
        iid_index = self._iid_range[iid_index_or_none]
        sid_index = self._sid_range[sid_index_or_none]

        iid_count_out = len(iid_index)
        sid_count_out = len(sid_index)

        if not force_python_only:
            from sgkit_plink import wrap_plink_parser

            val = np.zeros((iid_count_out, sid_count_out), order=order, dtype=dtype)
            bed_fn = open_bed._name_of_other_file(self.filename, "bed", "bed")

            num_threads = self._get_num_threads()

            if iid_count_in > 0 and sid_count_in > 0:
                if dtype == np.int8:
                    if order == "F":
                        wrap_plink_parser.readPlinkBedFile2int8FAAA(
                            bed_fn.encode("ascii"),
                            iid_count_in,
                            sid_count_in,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    elif order == "C":
                        wrap_plink_parser.readPlinkBedFile2int8CAAA(
                            bed_fn.encode("ascii"),
                            iid_count_in,
                            sid_count_in,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    else:
                        raise Exception(
                            "order '{0}' not known, only 'F' and 'C'".format(order)
                        )
                elif dtype == np.float64:
                    if order == "F":
                        wrap_plink_parser.readPlinkBedFile2doubleFAAA(  #!!!cmk double check that these check the format.If they don't, be sure checkformat is called sometime
                            bed_fn.encode("ascii"),
                            iid_count_in,
                            sid_count_in,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    elif order == "C":
                        wrap_plink_parser.readPlinkBedFile2doubleCAAA(
                            bed_fn.encode("ascii"),
                            iid_count_in,
                            sid_count_in,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    else:
                        raise Exception(
                            "order '{0}' not known, only 'F' and 'C'".format(
                                order
                            )  #!!!cmk or A
                        )
                elif dtype == np.float32:
                    if order == "F":
                        wrap_plink_parser.readPlinkBedFile2floatFAAA(
                            bed_fn.encode("ascii"),
                            iid_count_in,
                            sid_count_in,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    elif order == "C":
                        wrap_plink_parser.readPlinkBedFile2floatCAAA(
                            bed_fn.encode("ascii"),
                            iid_count_in,
                            sid_count_in,
                            self.count_A1,
                            iid_index,
                            sid_index,
                            val,
                            num_threads,
                        )
                    else:
                        raise Exception(
                            "order '{0}' not known, only 'F' and 'C'".format(order)
                        )
                else:
                    raise Exception(
                        "dtype '{0}' not known, only float64 and float32 (and sometimes int8)".format(
                            dtype
                        )
                    )

        else:
            if not self.count_A1:
                byteZero = 0
                byteThree = 2
            else:
                byteZero = 2
                byteThree = 0
            if dtype == np.int8:
                missing = -127
            else:
                missing = np.nan
            # An earlier version of this code had a way to read consecutive SNPs of code in one read. May want #!!!cmk understand these messages
            # to add that ability back to the code.
            # Also, note that reading with python will often result in non-contiguous memory, so the python standardizers will automatically be used, too.
            # logging.warn("using pure python plink parser (might be much slower!!)")
            val = np.zeros(
                ((int(np.ceil(0.25 * iid_count_in)) * 4), sid_count_out),
                order=order,
                dtype=dtype,
            )  # allocate it a little big

            bedfile = self._name_of_other_file(self.filename, "bed", "bed")
            with open(bedfile, "rb") as filepointer:
                for SNPsIndex, bimIndex in enumerate(sid_index):

                    startbit = int(np.ceil(0.25 * iid_count_in) * bimIndex + 3)
                    filepointer.seek(startbit)
                    nbyte = int(np.ceil(0.25 * iid_count_in))
                    bytes = np.array(bytearray(filepointer.read(nbyte))).reshape(
                        (int(np.ceil(0.25 * iid_count_in)), 1), order="F"
                    )

                    val[3::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[3::4, SNPsIndex : SNPsIndex + 1][bytes >= 64] = missing
                    val[3::4, SNPsIndex : SNPsIndex + 1][bytes >= 128] = 1
                    val[3::4, SNPsIndex : SNPsIndex + 1][bytes >= 192] = byteThree
                    bytes = np.mod(bytes, 64)
                    val[2::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[2::4, SNPsIndex : SNPsIndex + 1][bytes >= 16] = missing
                    val[2::4, SNPsIndex : SNPsIndex + 1][bytes >= 32] = 1
                    val[2::4, SNPsIndex : SNPsIndex + 1][bytes >= 48] = byteThree
                    bytes = np.mod(bytes, 16)
                    val[1::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[1::4, SNPsIndex : SNPsIndex + 1][bytes >= 4] = missing
                    val[1::4, SNPsIndex : SNPsIndex + 1][bytes >= 8] = 1
                    val[1::4, SNPsIndex : SNPsIndex + 1][bytes >= 12] = byteThree
                    bytes = np.mod(bytes, 4)
                    val[0::4, SNPsIndex : SNPsIndex + 1] = byteZero
                    val[0::4, SNPsIndex : SNPsIndex + 1][bytes >= 1] = missing
                    val[0::4, SNPsIndex : SNPsIndex + 1][bytes >= 2] = 1
                    val[0::4, SNPsIndex : SNPsIndex + 1][bytes >= 3] = byteThree
                val = val[iid_index, :]  # reorder or trim any extra allocation
                if not open_bed._array_properties_are_ok(val, order, dtype):
                    val = val.copy(order=order)

        return val

    @staticmethod
    def _split_index(index):
        if not isinstance(index, tuple):
            index = (None, index)
        iid_index = open_bed._fix_up_index(index[0])
        sid_index = open_bed._fix_up_index(index[1])
        return iid_index, sid_index

    @staticmethod
    def _fix_up_index(index):
        if index is None:  # make a shortcut for None
            return slice(None)
        try:  # If index is an int, return it in an array
            index = index.__index__()  # (see
            # https://stackoverflow.com/questions/3501382/checking-whether-a-variable-is-an-integer-or-not)
            return [index]
        except Exception:
            pass
        return index

    #!!!cmk kill
    # @staticmethod
    # def _fixup_input(
    #    input,
    #    count=None,
    #    empty_creator=_default_empty_creator,
    #    dtype=None,
    #    none_num_ok=False,
    # ):
    #    if none_num_ok and input is None:
    #        return input

    #    if input is None or len(input) == 0:
    #        input = empty_creator(count)
    #    elif not isinstance(input, np.ndarray):
    #        input = np.array(input, dtype=dtype)

    #    assert (
    #        count is None or len(input) == count
    #    ), "Expect length of {0} for input {1}".format(count, input)

    #    return input

    #!!!cmk put the methods in a good order

    @staticmethod
    def _write_fam_bim(basefilename, metadata, suffix_of_interest):
        filename = open_bed._name_of_other_file(basefilename, "bed", suffix_of_interest)

        fam_bim_list = []
        for (key, (suffix, position, _, _, _),) in open_bed._meta_meta:
            if suffix == suffix_of_interest:
                assert len(fam_bim_list) == position, "real assert"
                fam_bim_list.append(metadata[key])

        if suffix_of_interest == "fam":
            sep = " "
        else:
            sep = "\t"

        with open(filename, "w") as filepointer:
            for index in range(len(fam_bim_list[0])):
                filepointer.write(sep.join(seq[index] for seq in fam_bim_list) + "\n")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import os

    if True:
        import numpy as np
        from sgkit_plink._open_bed import open_bed

        # Can get file from https://www.dropbox.com/sh/xluk9opjiaobteg/AABgEggLk0ZoO0KQq0I4CaTJa?dl=0
        bigfile = r"M:\deldir\genbgen\2\merged_487400x220000.1.bed"
        # bigfile = '/mnt/m/deldir/genbgen/2/merged_487400x220000.1.bed'
        with open_bed(bigfile, num_threads=20) as bed:
            sid_batch = 22 * 1000
            for sid_start in range(0, 10 * sid_batch, sid_batch):
                slicer = np.s_[:10000, sid_start : sid_start + sid_batch]
                print(slicer)
                val = bed.read(slicer)
            print(val.shape)

    if False:
        file = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests\data/plink_sim_10s_100v_10pmiss.bed"  #!!!cmk remove absolute reference
        with open_bed(file) as bed:
            print(bed.iid)
            print(bed.shape)
            val = bed.read()
            print(val)

    if False:

        # bed_file = example_file('doc/ipynb/all.*','*.bed')
        bed_file = r"F:\backup\carlk4d\data\carlk\cachebio\genetics\onemil\id1000000.sid_1000000.seed0.byiid\iid990000to1000000.bed"
        bed = Bed(bed_file, count_A1=False)
        snpdata1 = bed[:, :1000].read()
        snpdata2 = bed[:, :1000].read(dtype="int8", _require_float32_64=False)
        print(snpdata2)
        snpdata3 = bed[:, :1000].read(
            dtype="int8", order="C", _require_float32_64=False
        )
        print(snpdata3)
        snpdata3.val = snpdata3.val.astype("float32")
        snpdata3.val.dtype

    if False:
        from sgkit_plink import Bed, SnpGen

        iid_count = 487409
        sid_count = 5000
        sid_count_max = 5765294
        sid_batch_size = 50

        sid_batch_count = -(sid_count // -sid_batch_size)
        sid_batch_count_max = -(sid_count_max // -sid_batch_size)
        snpgen = SnpGen(seed=234, iid_count=iid_count, sid_count=sid_count_max)

        for batch_index in range(sid_batch_count):
            sid_index_start = batch_index * sid_batch_size
            sid_index_end = (batch_index + 1) * sid_batch_size  # what about rounding
            filename = r"d:\deldir\rand\fakeukC{0}x{1}-{2}.bed".format(
                iid_count, sid_index_start, sid_index_end
            )
            if not os.path.exists(filename):  #!!!cmk use pathlib?
                Bed.write(
                    filename + ".temp", snpgen[:, sid_index_start:sid_index_end].read()
                )
                os.rename(filename + ".temp", filename)

    if False:
        from sgkit_plink import Pheno, Bed

        filename = r"m:\deldir\New folder (4)\all_chr.maf0.001.N300.bed"
        iid_count = 300
        iid = [["0", "iid_{0}".format(iid_index)] for iid_index in range(iid_count)]
        bed = Bed(filename, iid=iid, count_A1=False)
        print(bed.iid_count)

    if False:
        from pysnptools.util import example_file

        pheno_fn = example_file("pysnptools/examples/toydata.phe")

    if False:
        from sgkit_plink import Pheno, Bed

        print(os.getcwd())
        snpdata = Pheno("../examples/toydata.phe").read()  # Read data from Pheno format
        # pstutil.create_directory_if_necessary("tempdir/toydata.5chrom.bed")
        Bed.write(
            "tempdir/toydata.5chrom.bed", snpdata, count_A1=False
        )  # Write data in Bed format

    import doctest

    #!!!cmk put this back
    # doctest.testmod(
    #    optionflags=doctest.ELLIPSIS
    # )  #!!!cmk how do you doctest with PyTest?
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
