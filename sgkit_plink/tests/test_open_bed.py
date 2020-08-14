import numpy as np
from pathlib import Path
import pytest
from sgkit_plink._open_bed import open_bed
import logging  #!!!cmk how does sgkit do logging messages?


def test_read1():
    file = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests\data/plink_sim_10s_100v_10pmiss.bed"  #!!!cmk remove absolute reference
    with open_bed(file) as bed:
        assert bed.iid_count == 10
        assert bed.fid[-1] == "0"
        assert bed.iid[-1] == "9"
        assert bed.shape == (10, 100)
        val = bed.read()
        assert (
            val.mean() == -13.142
        )  # really shouldn't do mean on data where -127 represents missing
        assert bed.chromosome[-1] == "1"
        assert bed.bp_position[-1] == 100
        #!!!cmk test reading into other dtypes


#!!!cmk write test
# def test_write():
#    in_file = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests\data/plink_sim_10s_100v_10pmiss.bed"  #!!!cmk remove absolute reference
#    out_file = r"m:/deldir/out.bed"  #!!!cmk remove absolute reference
#    with open_bed(in_file) as bed:
#        val0 = bed.read()
#        old_iid = np.array([bed.fid, bed.iid]).T
#        pos = np.array(
#            [bed.chromosome.astype("int"), bed.cm_position, bed.bp_position]
#        ).T
#        open_bed.write(
#            out_file, val0, iid=old_iid, sid=bed.sid, pos=pos,
#        )
#        with open_bed(out_file) as bed1:
#            assert (val0 == bed1.read()).all()  #!!!cmk use array_equal
#            assert (bed.iid == bed1.iid).all()
#            assert (bed.sid == bed1.sid).all()
#            assert (
#                bed.chromosome.astype("float") == bed1.chromosome.astype("float")
#            ).all()  #!!!cmk remove the 'astype('float')'
#            assert (bed.cm_position == bed1.cm_position).all()
#            assert (bed.bp_position == bed1.bp_position).all()

#    val_float = val0.astype("float")
#    val_float[0, 0] = 0.5
#    for force_python_only in [
#        True
#    ]:  #!!!cmk It is a bug that the C++ version doesn't catch e.g. .5 as input
#        with pytest.raises(ValueError):
#            open_bed.write(
#                out_file,
#                val_float,
#                iid=old_iid,
#                sid=bed.sid,
#                pos=pos,
#                force_python_only=force_python_only,
#            )


def test_overrides():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    with open_bed(base / "data/distributed_bed_test1_X.bed") as bed:
        fid = bed.fid
        iid = bed.iid
        father = bed.father
        mother = bed.mother
        sex = bed.sex
        pheno = bed.pheno
        chromosome = bed.chromosome
        sid = bed.sid
        cm_position = bed.cm_position
        bp_position = bed.bp_position
        allele_1 = bed.allele_1
        allele_2 = bed.allele_2
    # lock in the expected results: np.savez(base / "data/distributed_bed_test1_X.metadata.npz",fid=fid,iid=iid,father=father,mother=mother,sex=sex,pheno=pheno,chromosome=chromosome,sid=sid,cm_position=cm_position,bp_position=bp_position,allele_1=allele_1,allele_2=allele_2)
    property_dict = np.load(base / "data/distributed_bed_test1_X.metadata.npz")
    assert np.array_equal(property_dict["fid"], fid)
    assert np.array_equal(property_dict["iid"], iid)
    assert np.array_equal(property_dict["father"], father)
    assert np.array_equal(property_dict["mother"], mother)
    assert np.array_equal(property_dict["sex"], sex)
    assert np.array_equal(property_dict["pheno"], pheno)
    assert np.array_equal(property_dict["chromosome"], chromosome)
    assert np.array_equal(property_dict["sid"], sid)
    assert np.array_equal(property_dict["cm_position"], cm_position)
    assert np.array_equal(property_dict["bp_position"], bp_position)
    assert np.array_equal(property_dict["allele_1"], allele_1)
    assert np.array_equal(property_dict["allele_2"], allele_2)

    with pytest.raises(KeyError):
        open_bed(
            base / "data/distributed_bed_test1_X.bed", metadata={"unknown": [3, 4, 4]}
        )
    with open_bed(
        base / "data/distributed_bed_test1_X.bed", metadata={"iid": None}
    ) as bed1:
        assert np.array_equal(bed1.iid, property_dict["iid"])
    with open_bed(
        base / "data/distributed_bed_test1_X.bed", metadata={"iid": []}
    ) as bed1:
        assert bed1.iid.dtype.type is np.str_
        assert len(bed1.iid) == 0
        with pytest.raises(ValueError):
            bed1.father

    with open_bed(
        base / "data/distributed_bed_test1_X.bed",
        metadata={"sid": [i for i in range(len(sid))]},
    ) as bed1:
        assert bed1.sid.dtype.type is np.str_
        assert bed1.sid[0] == "0"
    with pytest.raises(ValueError):
        open_bed(
            base / "data/distributed_bed_test1_X.bed",
            metadata={"sex": ["F" for i in range(len(sex))]},
        )  # Sex must be coded as a number
    with open_bed(
        base / "data/distributed_bed_test1_X.bed",
        metadata={"sid": np.array([i for i in range(len(sid))])},
    ) as bed1:
        assert bed1.sid.dtype.type is np.str_
        assert bed1.sid[0] == "0"
    with pytest.raises(ValueError):
        open_bed(
            base / "data/distributed_bed_test1_X.bed",
            metadata={"sid": np.array([(i, i) for i in range(len(sid))])},
        )
    with open_bed(
        base / "data/distributed_bed_test1_X.bed", metadata={"sid": [1, 2, 3]}
    ) as bed1:
        with pytest.raises(ValueError):
            bed1.chromosome


def test_str():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    with open_bed(base / "data/distributed_bed_test1_X.bed") as bed:
        assert "open_bed(" in str(bed)


def test_bad_bed():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    with pytest.raises(ValueError):
        open_bed(base / "data/badfile.bed")


# def test_read_empty_metafiles():
#    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
#    #!!!need to generate empty bed file with zero length fam and bim files for testing and then test it
#    with open_bed(base / "data/distributed_bed_test1_X.bed") as bed:


#!!!cmk rather slow
def test_properties():
    file = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests\data/plink_sim_10s_100v_10pmiss.bed"  #!!!cmk remove absolute reference
    with open_bed(file) as bed:
        iid_list = bed.iid.tolist()
        sid_list = bed.sid.tolist()
        chromosome_list = bed.chromosome.tolist()
    test_count = 0
    for iid in [None, iid_list, np.array(iid_list)]:
        for iid_count in [None, len(iid_list)]:
            for iid_before_read in [False, True]:
                for iid_after_read in [False, True]:
                    for sid in [None, sid_list, np.array(sid_list)]:
                        for sid_count in [None, len(sid_list)]:
                            for sid_before_read in [False, True]:
                                for sid_after_read in [False, True]:
                                    for chromosome in [
                                        None,
                                        chromosome_list,
                                        np.array(chromosome_list),
                                    ]:
                                        for chromosome_before_read in [False, True]:
                                            for chromosome_after_read in [False, True]:
                                                with open_bed(
                                                    file,
                                                    iid_count=iid_count,
                                                    sid_count=sid_count,
                                                    metadata={
                                                        "iid": iid,
                                                        "sid": sid,
                                                        "chromosome": chromosome,
                                                    },
                                                ) as bed:
                                                    logging.info(f"Test {test_count}")
                                                    test_count += 1
                                                    if iid_before_read:
                                                        assert np.array_equal(
                                                            bed.iid, iid_list
                                                        )
                                                    if sid_before_read:
                                                        assert np.array_equal(
                                                            bed.sid, sid_list
                                                        )
                                                    if chromosome_before_read:
                                                        assert np.array_equal(
                                                            bed.chromosome,
                                                            chromosome_list,
                                                        )
                                                    val = bed.read()
                                                    assert val.shape == (
                                                        len(iid_list),
                                                        len(sid_list),
                                                    )
                                                    if iid_after_read:
                                                        assert np.array_equal(
                                                            bed.iid, iid_list
                                                        )
                                                    if sid_after_read:
                                                        assert np.array_equal(
                                                            bed.sid, sid_list
                                                        )
                                                    if chromosome_after_read:
                                                        assert np.array_equal(
                                                            bed.chromosome,
                                                            chromosome_list,
                                                        )
                                                    # bed._assert_iid_sid_chromosome()


def test_c_reader_bed():
    for force_python_only in [False, True]:
        base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
        bed = open_bed(
            base / "data/distributed_bed_test1_X.bed", count_A1=False
        )  # !!!cmk improve the name of this test file (it contains missing)

        val = bed.read(order="F", dtype="float64", force_python_only=force_python_only)
        assert val.dtype.type == np.float64
        ref_val = reference_val()
        ref_val = ref_val * -1 + 2
        assert np.allclose(ref_val, val, rtol=1e-05, atol=1e-05, equal_nan=True)

        val = bed.read(order="F", force_python_only=False)
        assert val.dtype.type == np.int8
        ref_val[ref_val != ref_val] = -127
        ref_val = ref_val.astype("int8")
        ref_val = ref_val.astype("int8")
        assert np.all(ref_val == val)

        bed.close()

        base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
        with open_bed(base / "data/distributed_bed_test1_X.bed") as bed:
            val = bed.read(
                order="F", dtype="float64", force_python_only=force_python_only
            )
            ref_val = reference_val()
            assert np.allclose(ref_val, val, rtol=1e-05, atol=1e-05, equal_nan=True)


def reference_val():  #!!!cmk fix this so not loading over and over again
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    val = np.load(base / "data/distributed_bed_test1_X.val.npy")
    return val


def test_bed_int8():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    with open_bed(base / "data/distributed_bed_test1_X.bed") as bed:
        for force_python_only in [False, True]:
            for order in ["F", "C"]:
                val = bed.read(
                    dtype="int8", force_python_only=force_python_only, order=order
                )
                assert val.dtype == "int8"
                assert (val.flags["C_CONTIGUOUS"] and order == "C") or (
                    val.flags["F_CONTIGUOUS"] and order == "F"
                )
                ref_val = reference_val()
                ref_val[ref_val != ref_val] = -127
                ref_val = ref_val.astype("int8")
                assert np.array_equal(ref_val, val)
                #!!!cmk add write test later
                # output = "tempdir/snpreader/int8.bed" #!!!cmk where does output go?
                # create_directory_if_necessary(output)
                # for count_A1 in [False,True]:
                #    bed2 = Bed.write(output,snpdata,count_A1=count_A1,_require_float32_64=False,force_python_only=force_python_only)
                #    assert np.allclose(bed2.read(dtype='int8',_require_float32_64=False,force_python_only=force_python_only).val, ref.val, equal_nan=True)


# !!!cmk put this test back in?
# def test_too_slow_write_bedbig():
#    iid_count = 100000
#    sid_count = 50000
#    from pysnptools.snpreader import SnpData
#    iid = np.array([[str(i),str(i)] for i in range(iid_count)])
#    sid = np.array(["sid_{0}".format(i) for i in range(sid_count)])
#    pos = np.array([[i,i,i] for i in range(sid_count)])
#    np.random.seed(0)
#    snpdata = SnpData(iid,sid,np.zeros((iid_count,sid_count)),pos=pos) #random.choice((0.0,1.0,2.0,float("nan")),size=(iid_count,sid_count)))
#    output = "tempdir/bedbig.{0}.{1}".format(iid_count,sid_count)
#    create_directory_if_necessary(output)
#    Bed.write(output, snpdata, count_A1=False)
#    snpdata2 = Bed(output,count_A1=False).read()
#    np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)


def test_write_bed_f64cpp():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    with open_bed(base / "data/distributed_bed_test1_X.bed") as bed:
        for iid_index in [0, 1, 5]:
            for force_python_only in [False, True]:
                val = bed.read(
                    np.s_[0:iid_index, :],
                    order="F",
                    dtype=np.float64,
                    force_python_only=force_python_only,
                )
                assert val.shape == (iid_index, 100)
                #!!!cmk add write test later
                # output = "tempdir/toydata.F64cpp.{0}".format(iid_index)
                # create_directory_if_necessary(output)
                # Bed.write(output, snpdata ,count_A1=False)
                # snpdata2 = Bed(output,count_A1=False).read()
                # np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)


#!!!cmk write test
# def test_write_x_x_cpp():
#    base = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests"
#    for count_A1 in [False, True]:
#        snpreader = Bed(self.currentFolder + "/examples/toydata.5chrom.bed",count_A1=count_A1)
#        for order in ['C','F']:
#            for dtype in [np.float32,np.float64]:
#                snpdata = snpreader.read(order=order,dtype=dtype)
#                snpdata.val[-1,0] = float("NAN")
#                output = "tempdir/toydata.{0}{1}.cpp".format(order,"32" if dtype==np.float32 else "64")
#                create_directory_if_necessary(output)
#                Bed.write(output, snpdata, count_A1=count_A1)
#                snpdata2 = Bed(output,count_A1=count_A1).read()
#                np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)


def test_respect_read_inputs():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    ref_val_float = reference_val()

    ref_val_int8 = ref_val_float.astype("int8")
    ref_val_int8[ref_val_float != ref_val_float] = -127

    with open_bed(base / "data/distributed_bed_test1_X.bed") as bed:
        for order in ["F", "C", "A"]:
            for dtype in [np.int8, np.float32, np.float64]:
                for force_python_only in [True, False]:
                    val = bed.read(
                        order=order, dtype=dtype, force_python_only=force_python_only
                    )
                    has_right_order = (
                        order == "A"
                        or (order == "C" and val.flags["C_CONTIGUOUS"])
                        or (order == "F" and val.flags["F_CONTIGUOUS"])
                    )
                    assert val.dtype == dtype and has_right_order
                    ref_val = ref_val_int8 if dtype == np.int8 else ref_val_float
                    assert np.allclose(ref_val, val, equal_nan=True)


def test_threads():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    ref_val_float = reference_val()

    ref_val_float = reference_val()
    ref_val_int8 = ref_val_float.astype("int8")
    ref_val_int8[ref_val_float != ref_val_float] = -127

    for num_threads in [1, 4]:
        with open_bed(
            base / "data/distributed_bed_test1_X.bed", num_threads=num_threads
        ) as bed:
            val = bed.read()
            assert np.allclose(ref_val_int8, val, equal_nan=True)


#!!!cmk add write tests
# def test_writes():
#    base = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests"
#    from pysnptools.snpreader import SnpData, SnpHdf5, SnpNpz, SnpMemMap

#    the_class_and_suffix_list = [
#                                    (Bed,"bed",lambda filename:Bed(filename,count_A1=False),None),
#                                ]
#    cant_do_col_prop_none_set = {"dense","distributed_bed"}
#    cant_do_col_len_0_set = {"distributed_bed"}
#    cant_do_row_count_zero_set = {'dense','ped','pheno'}
#    can_swap_0_2_set = {'ped'}
#    can_change_col_names_set = {'pheno'}
#    ignore_fam_id_set = {'dense'}
#    ignore_pos_set = {'dense','pheno'}
#    erase_any_write_dir = {'distributed_bed'}


#    #===================================
#    #    Starting main function
#    #===================================
#    logging.info("starting 'test_writes'")
#    np.random.seed(0)
#    output_template = "tempdir/snpreader/writes.{0}.{1}"
#    create_directory_if_necessary(output_template.format(0,"npz"))
#    i = 0
#    for row_count in [0,5,2,1]:
#        for col_count in [4,2,1,0]:
#            val = np.random.randint(0,4,size=(row_count,col_count))*1.0
#            val[val==3]=np.NaN
#            row = [('0','0'),('1','1'),('2','2'),('3','3'),('4','4')][:row_count]
#            col = ['s0','s1','s2','s3','s4'][:col_count]
#            for is_none in [True,False]:
#                row_prop = None
#                col_prop = None if is_none else [(x,x,x) for x in range(5)][:col_count]
#                snpdata = SnpData(iid=row,sid=col,val=val,pos=col_prop,name=str(i))
#                for the_class,suffix,constructor,writer in the_class_and_suffix_list:
#                    constructor = constructor or (lambda filename: the_class(filename))
#                    writer = writer or (lambda filename,_data: the_class.write(filename,_data))

#                    if col_count == 0 and suffix in cant_do_col_len_0_set:
#                        continue
#                    if col_prop is None and suffix in cant_do_col_prop_none_set:
#                        continue
#                    if row_count==0 and suffix in cant_do_row_count_zero_set:
#                        continue
#                    filename = output_template.format(i,suffix)
#                    logging.info(filename)
#                    i += 1
#                    if suffix in erase_any_write_dir and os.path.exists(filename):
#                        shutil.rmtree(filename)
#                    ret = writer(filename,snpdata)
#                    assert ret is not None
#                    for subsetter in [None, np.s_[::2,::3]]:
#                        reader = constructor(filename)
#                        _fortesting_JustCheckExists().input(reader)
#                        subreader = reader if subsetter is None else reader[subsetter[0],subsetter[1]]
#                        readdata = subreader.read(order='C')
#                        expected = snpdata if subsetter is None else snpdata[subsetter[0],subsetter[1]].read()
#                        if not suffix in can_swap_0_2_set:
#                            assert np.allclose(readdata.val,expected.val,equal_nan=True)
#                        else:
#                            for col_index in range(readdata.col_count):
#                                assert (np.allclose(readdata.val[:,col_index],expected.val[:,col_index],equal_nan=True) or
#                                        np.allclose(readdata.val[:,col_index]*-1+2,expected.val[:,col_index],equal_nan=True))
#                        if not suffix in ignore_fam_id_set:
#                            assert np.array_equal(readdata.row,expected.row)
#                        else:
#                            assert np.array_equal(readdata.row[:,1],expected.row[:,1])
#                        if not suffix in can_change_col_names_set:
#                            assert np.array_equal(readdata.col,expected.col)
#                        else:
#                            assert readdata.col_count==expected.col_count
#                        assert np.array_equal(readdata.row_property,expected.row_property) or (readdata.row_property.shape[1]==0 and expected.row_property.shape[1]==0)

#                        if not suffix in ignore_pos_set:
#                            assert np.allclose(readdata.col_property,expected.col_property,equal_nan=True) or (readdata.col_property.shape[1]==0 and expected.col_property.shape[1]==0)
#                        else:
#                            assert len(readdata.col_property)==len(expected.col_property)
#                    try:
#                        os.remove(filename)
#                    except:
#                        pass
#    logging.info("done with 'test_writes'")


def test_index():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    ref_val_float = reference_val()

    ref_val_float = reference_val()
    ref_val_int8 = ref_val_float.astype("int8")
    ref_val_int8[ref_val_float != ref_val_float] = -127

    with open_bed(
        base / "data/distributed_bed_test1_X.bed"
    ) as bed:  #!!!cmk maybe repeat with a 2nd file that doesn't have so many 2's in it and isn't square
        val = bed.read()
        assert np.allclose(ref_val_int8, val, equal_nan=True)

        val = bed.read(2)
        assert np.allclose(ref_val_int8[:, [2]], val, equal_nan=True)

        val = bed.read((2))
        assert np.allclose(ref_val_int8[:, [2]], val, equal_nan=True)

        val = bed.read((None, 2))
        assert np.allclose(ref_val_int8[:, [2]], val, equal_nan=True)

        val = bed.read((1, 2))
        assert np.allclose(ref_val_int8[[1], [2]], val, equal_nan=True)

        val = bed.read([2, -2])
        assert np.allclose(ref_val_int8[:, [2, -2]], val, equal_nan=True)

        val = bed.read(([1, -1], [2, -2]))
        assert np.allclose(ref_val_int8[[1, -1], :][:, [2, -2]], val, equal_nan=True)

        iid_bool = ([False, False, True] * bed.iid_count)[: bed.iid_count]
        sid_bool = ([True, False, True] * bed.sid_count)[: bed.sid_count]
        val = bed.read(sid_bool)
        assert np.allclose(ref_val_int8[:, sid_bool], val, equal_nan=True)

        val = bed.read((iid_bool, sid_bool))
        assert np.allclose(ref_val_int8[iid_bool, :][:, sid_bool], val, equal_nan=True)

        val = bed.read((1, sid_bool))
        assert np.allclose(ref_val_int8[[1], :][:, sid_bool], val, equal_nan=True)

        slicer = np.s_[::2, ::3]
        val = bed.read(slicer[1])
        assert np.allclose(ref_val_int8[:, slicer[1]], val, equal_nan=True)

        val = bed.read(slicer)
        assert np.allclose(ref_val_int8[slicer], val, equal_nan=True)

        val = bed.read((1, slicer[1]))
        assert np.allclose(ref_val_int8[[1], slicer[1]], val, equal_nan=True)


def test_shape():
    base = Path(r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests")
    with open_bed(base / "data/plink_sim_10s_100v_10pmiss.bed") as bed:
        assert bed.shape == (10, 100)


if __name__ == "__main__":  #!!cmk is this wanted?
    logging.basicConfig(level=logging.INFO)

    test_read1()
    pytest.main([__file__])
