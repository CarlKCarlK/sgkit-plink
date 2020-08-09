import numpy as np
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
        assert val.mean() == -13.142 # really shouldn't do mean on data where -127 represents missing
        assert bed.chromosome[-1] == "1"
        assert bed.bp_position[-1] == 100
        #!!!cmk test reading into other dtypes

def test_write():
    in_file = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests\data/plink_sim_10s_100v_10pmiss.bed"  #!!!cmk remove absolute reference
    out_file = r"m:/deldir/out.bed"  #!!!cmk remove absolute reference
    with open_bed(in_file) as bed:
        val0 = bed.read()
        old_iid = np.array(
            [bed.fid, bed.iid]
        ).T
        pos = np.array(
            [bed.chromosome.astype("int"), bed.cm_position, bed.bp_position]
        ).T
        open_bed.write(
            out_file, val0, iid=old_iid, sid=bed.sid, pos=pos,
        )
        with open_bed(out_file) as bed1:
            assert (val0 == bed1.read()).all()  #!!!cmk use array_equal
            assert (bed.iid == bed1.iid).all()
            assert (bed.sid == bed1.sid).all()
            assert (
                bed.chromosome.astype("float") == bed1.chromosome.astype("float")
            ).all()  #!!!cmk remove the 'astype('float')'
            assert (bed.cm_position == bed1.cm_position).all()
            assert (bed.bp_position == bed1.bp_position).all()

    val_float = val0.astype("float")
    val_float[0, 0] = 0.5
    for force_python_only in [True]: #!!!cmk It is a bug that the C++ version doesn't catch e.g. .5 as input
        with pytest.raises(ValueError):
            open_bed.write(
                out_file,
                val_float,
                iid=old_iid,
                sid=bed.sid,
                pos=pos,
                force_python_only=force_python_only,
            )


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
                                                    iid=iid,
                                                    sid=sid,
                                                    chromosome=chromosome,
                                                ) as bed:
                                                    logging.info(f"Test {test_count}")
                                                    test_count += 1
                                                    if iid_before_read:
                                                        assert np.array_equal(bed.iid, iid_list)
                                                    if sid_before_read:
                                                        assert np.array_equal(bed.sid, sid_list)
                                                    if chromosome_before_read:
                                                        assert np.array_equal(
                                                            bed.chromosome, chromosome_list
                                                        )
                                                    val = bed.read()
                                                    assert val.shape == (
                                                        len(iid_list),
                                                        len(sid_list),
                                                    )
                                                    if iid_after_read:
                                                        assert np.array_equal(bed.iid, iid_list)
                                                    if sid_after_read:
                                                        assert np.array_equal(bed.sid, sid_list)
                                                    if chromosome_after_read:
                                                        assert np.array_equal(
                                                            bed.chromosome, chromosome_list
                                                        )
                                                    # bed._assert_iid_sid_chromosome()


def test_shape():
    with open_bed(r'D:\OneDrive\programs\pstsgkit\tests\datasets\all_chr.maf0.001.N300.bed') as bed:
        assert bed.shape==(300,1015)


if __name__ == "__main__":  #!!cmk is this wanted?
    logging.basicConfig(level=logging.INFO)

    test_properties()  #!!!cmk
    pytest.main([__file__])
