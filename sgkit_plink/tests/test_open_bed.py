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
        assert bed.pos

def test_write():



if __name__ == "__main__": #!!cmk is this wanted?
    pytest.main([__file__])
