import pytest

from app.utils.helpers import validate_chromosome


def test_validate_chromosome():
    assert validate_chromosome(chr=1, prefix=None)
    assert validate_chromosome(chr="chr1", prefix="chr")
    assert validate_chromosome(chr=23, prefix=None, x_y_numeric=True)
    assert validate_chromosome(chr="chX", prefix="ch", x_y_numeric=False)

    with pytest.raises(ValueError):
        validate_chromosome(chr=1, prefix="ch")
    with pytest.raises(ValueError):
        validate_chromosome(chr="chr1", prefix="None")
    with pytest.raises(ValueError):
        validate_chromosome(chr=23, prefix=None, x_y_numeric=False)
    with pytest.raises(ValueError):
        validate_chromosome(chr="chX", prefix="ch", x_y_numeric=True)
