from typing import Any, List


def validate_chromosome(
    chr: str | int, prefix: str | None = "chr", x_y_numeric: bool = False
) -> bool:
    """Validate that the chromosome is in the correct format

    :param chr: The chromosome representation, Examples: chr1, 1, chrX, chrome23
    :type chr: str | int
    :param prefix: The part of `chr` that is a string, defaults to "chr"
    :type prefix: str | None, optional
    :param x_y_numeric: Whether X and Y are represented as a number, defaults to False
    :type x_y_numeric: bool, optional
    :raises ValueError: If the validation fails
    :return: The validation status (True if passing, otherwise a ValueError is raised)
    :rtype: bool
    """
    chrs: List[str] = [str(c) for c in range(1, 23)]
    chr = str(chr)
    if x_y_numeric is True:
        chrs.append("23")
    else:
        chrs.extend(["X", "Y"])
    if prefix is not None:
        if not chr.startswith(prefix):
            raise ValueError(f"Chromosome must start with {prefix}.")
        chr = chr.replace(prefix, "")
    if chr not in chrs:
        raise ValueError(f"Chromosome must be one of {', '.join(chrs)}")

    return True
