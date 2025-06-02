def db2pow(db_val: float) -> float:
    """Convert dB value to power (linear scale)"""
    return 10**(db_val / 10)
