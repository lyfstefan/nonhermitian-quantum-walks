"""2026 reimplementation of PRA 102, 022218 (2020), Eqs. (1)-(2)."""
from .model import Walk, bulk_discriminant, bulk_regime

__all__ = ["Walk", "bulk_discriminant", "bulk_regime"]
