import pytest
from sgkit_plink._open_bed import open_bed

def test_read1():
    file = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests\data/plink_sim_10s_100v_10pmiss.bed"  #!!!cmk remove absolute reference
    with open_bed(file) as bed:
        assert bed.iid_count == 10
        assert (bed.iid[-1,:] == ['0','9']).all()
        assert bed.shape == (10,100)
        val = bed.read(force_python_only=True)
        assert val.mean()==-13.142
        assert (bed.pos[-1,:] == [1,0,100]).all()
        #!!!cmk test reading into other dtypes

#!!!cmk show example of reading where chrom==5
#!!!cmk could subsetting for iid, etc be introduced without need to create snpdata like object? Does the Bed file offer fast access to metadata that would make this worthwhile?

def test_write():
    in_file = r"D:\OneDrive\programs\sgkit-plink\sgkit_plink\tests\data/plink_sim_10s_100v_10pmiss.bed"  #!!!cmk remove absolute reference
    out_file = r"m:/deldir/out.bed"   #!!!cmk remove absolute reference
    with open_bed(in_file) as bed:
        val0 = val=bed.read(force_python_only=True)
        open_bed.write(out_file,val0,iid=bed.iid,sid=bed.sid,pos=bed.pos,force_python_only=True)
        with open_bed(out_file) as bed1:
            assert (val0 == bed1.read()).all() #!!!cmk use array_equal
            assert (bed.iid == bed1.iid).all()
            assert (bed.sid == bed1.sid).all()
            assert (bed.pos == bed1.pos).all()

    val_float = val0.astype('float')
    val_float[0,0]=0.5
    with pytest.raises(ValueError):
        open_bed.write(out_file,val_float,iid=bed.iid,sid=bed.sid,pos=bed.pos,force_python_only=True) #!!!cmk test on force_python=False, too


if __name__ == "__main__": #!!cmk is this wanted?
    pytest.main([__file__])
